import sys
import argparse
import os
import subprocess
import re
from collections import defaultdict
from multiprocessing.dummy import Pool as ThreadPool

SAM_POS_COL_INDEX = 3
SAM_CIGAR_COL_INDEX = 5


def get_bam_files_in_folder(bam_folder):
    """Get a list of all the bam file paths in a given folder"""
    bam_file_paths = []
    for root, dirs, files in os.walk(bam_folder):
        for file in files:
            if file.endswith(".bam"):
                 bam_file_paths.append(os.path.join(root, file))
    return bam_file_paths


def length_of_first_intronic_section(cigar_string):
    """Calculate the intron span's length.
    Intron length will be  last regex integer matching string before 'N' character. Example CIGAR string: '3M1D40M20N'
    """
    return int(re.findall('\d+', cigar_string.split('N')[0])[-1])


def num_of_matches_before_first_intronic_section(cigar_string):
    """Calculate the number of matches before intronic section of the CIGAR string. Example CIGAR string: '3M1D40M20N'
    """
    pre_N = cigar_string.split('N')[0]
    numbers = [int(regex_num) for regex_num in re.findall('\d+', pre_N)][0:-1]
    letters = [regex_letter for regex_letter in re.findall('[A-Z]', pre_N)]
    sum = 0
    for t in zip(numbers, letters):
        if t[1] not in {'I', 'S', 'P', 'H'}:
            sum += t[0]
    return sum


def write_cigar_debugging_info(cigar_string, num_matches_before_start, intron_length, pos, intron_start, intron_end):
    sys.stdout.write('>> {} ({} match, {} intronic)\t'
                     'read_start: {}\tintron_start: {}\tintron_end: {}\n'
                     .format(cigar_string,
                             num_matches_before_start,
                             intron_length,
                             pos,
                             intron_start,
                             intron_end))


def get_first_splice_junction(cigar_string, t_chrom, pos, verbose=False):
    intron_length = length_of_first_intronic_section(cigar_string)
    num_matches_before_intron_start = num_of_matches_before_first_intronic_section(cigar_string)
    intron_start_position = pos + num_matches_before_intron_start # inclusive
    intron_end_position = pos + num_matches_before_intron_start + intron_length # not inclusive

    if verbose:
        # For debugging purposes, if verbose is true, output some intermediate information
        write_cigar_debugging_info(cigar_string, num_matches_before_intron_start, intron_length,
                                   pos, intron_start_position, intron_end_position)

    unique_splice_junction = '{},{},{}'.format(t_chrom, intron_start_position, intron_end_position)
    return unique_splice_junction


def get_second_splice_junction(cigar_string, t_chrom, pos, first_splice_junction_end, verbose=False):
    section_after_first_intron = cigar_string.split('N')[1]
    num_matches_before_second_intron = \
        num_of_matches_before_first_intronic_section(section_after_first_intron)
    second_intron_length = length_of_first_intronic_section(section_after_first_intron)
    second_intron_start_position = first_splice_junction_end + num_matches_before_second_intron
    second_intron_end_position = second_intron_start_position + second_intron_length

    # For debugging purposes, if verbose is true, output some intermediate information
    if verbose:
        write_cigar_debugging_info(cigar_string, num_matches_before_second_intron,
                                   second_intron_length, pos, second_intron_start_position,
                                   second_intron_end_position)

    unique_splice_junction = '{},{},{}'.format(t_chrom, second_intron_start_position,
                                               second_intron_end_position)
    return unique_splice_junction


def get_end_position_from_junction(unique_splice_junction):
    return int(unique_splice_junction.split(',')[2])


def find_splice_junctions(bam_file_path, t_chrom, t_start, t_stop, verbose=False):
    """The meat of the script. Use samtools view along with some included filtering parameters described below to
    view reads and alignments representing intronic regions, then parse the CIGAR strings to figure out the start and
    end positions of these regions.

    # Samtools view options used here:
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # -F INT: Do not output alignments with any bits set in INT present in the FLAG field
    # -q INT: Skip alignments with MAPQ smaller than INT
    """
    t_location = '{}:{}-{}'.format(t_chrom, t_start, t_stop)
    pos = 'chr{}:{}-{}'.format(t_chrom, t_start, t_stop)
    samtools_view = 'samtools view -F 256 -F 512 -F 1024 -F 2048 -q 60 {} {} {}'.format(bam_file_path, pos,  t_location)
    # Run the samtools view command using subprocess, and place the stdout output into the sam_view variable as a string
    sam_view = subprocess.check_output(samtools_view.split())

    # Split the sam_view output string into each line comprising it
    sam_lines = str(sam_view).split('\n')

    # A dictionary that will keep track of how many times a particular junction is found
    splice_junctions = defaultdict(int)

    # Only look at the exon junction spanning lines, for which the CIGAR string will contain an 'N' character
    for line in sam_lines:
        try:
            split_line = line.split('\t')
            cigar_string = split_line[SAM_CIGAR_COL_INDEX]
            pos = int(split_line[SAM_POS_COL_INDEX])
            if 'N' in cigar_string and t_start < pos < t_stop:
                    # Account for only one intronic region, or only the first of two in the case there are two
                    if cigar_string.count('N') <= 2:
                        unique_splice_junction = get_first_splice_junction(cigar_string, t_chrom, pos, verbose=verbose)
                        splice_junctions[unique_splice_junction] += 1

                    # Account for very small retained introns (they must be small for the whole intron to have been
                    # captured within the length of one read), e.g. 13M221400N34M2658N29M
                    if cigar_string.count('N') == 2 and unique_splice_junction:
                        intron_end_position = get_end_position_from_junction(unique_splice_junction)  # not inclusive
                        unique_splice_junction = get_second_splice_junction(cigar_string, t_chrom, pos,
                                                                            intron_end_position, verbose=verbose)

                        splice_junctions[unique_splice_junction] += 1

        except IndexError:
            pass

    return splice_junctions


def get_id_from_bam_name(bam_file_path):
    """Strip away the .bam suffix and return the sample ID"""
    filename, file_extension = os.path.splitext(os.path.basename(bam_file_path))
    return filename


def summarize_splice_junctions(gene, global_event_counts, sample_ids, output_dir='.'):
    """Given a dictionary summarizing the global and per sample occurrence counts of unique splice junctions,
    output a file summarizing the results"""
    output_file_name = '{}/{}_{}_junctions.txt'.format(output_dir, gene, len(global_event_counts))
    with open(output_file_name, 'w') as f:
        ordered_events = sorted(global_event_counts.keys())
        for event in ordered_events:
            per_sample_counts = global_event_counts[event]
            gene, gene_type, chrom, start, end = event.split(',')

            total_count = 0

            sample_counts = [str(per_sample_counts.get(s, 0)) for s in sample_ids]

            for count in sample_counts:
                total_count += int(count)

            num_samples_with_this_event = len(per_sample_counts)

            entry = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(gene, gene_type, chrom, start, end, total_count,
                                                               num_samples_with_this_event, '\t'.join(sample_counts))
            f.write(entry)


def find_splice_junctions_for_gene(pool_arguments):
    """Find the splice junctions across all samples for one gene, writing the results to a file
    specific to this gene"""
    pa = pool_arguments
    # A dictionary to keep track of, for each unique splice junction event, how many times it appears in each of the
    # samples being considered as part of this analysis for this gene.
    global_event_counts = defaultdict(lambda: defaultdict())
    sample_id_to_bam_file_path = pa.get('sample_id_to_bam_file_path')
    sorted_sample_ids = sorted(sample_id_to_bam_file_path.keys())

    for sample_id in sorted_sample_ids:
        bam_file_path = sample_id_to_bam_file_path.get(sample_id)
        if pa.get('verbose'):
            sys.stdout.write("Sample: {}\tGene: {}\n".format(sample_id, pa.get('t_gene')))
        splice_junctions = find_splice_junctions(bam_file_path,
                                                 pa.get('t_chrom'),
                                                 int(pa.get('t_start')),
                                                 int(pa.get('t_stop')),
                                                 verbose=pa.get('verbose'))
        for splice_junction, count in splice_junctions.items():
            splice_junction_with_gene_info = '{},{},{}'.format(pa.get('t_gene'), pa.get('t_gene_type'), splice_junction)
            global_event_counts[splice_junction_with_gene_info][sample_id] = count
    summarize_splice_junctions(pa.get('t_gene'), global_event_counts, sorted_sample_ids, pa.get('output_dir'))


def map_splice_junction_discovery_across_genes(threads, transcript_file, sample_id_to_bam_file_path, verbose=False, output_dir=None):
    """Set up the parameters to map the splice junction discovery functions across all genes in a threaded manner"""
    pool_args = []
    pool = ThreadPool(threads)
    transcript_lines = open(transcript_file).readlines()
    for t in transcript_lines:
        if verbose:
            sys.stdout.write('>> {}\n'.format(t))
        if 'CHROM' in t:
            continue
        t_gene, t_enst, _, t_chrom, t_start, t_stop, t_gene_type = t.split('\t')
        t_gene = t_gene.split('\n')[0]
        t_gene_type = t_gene_type.split('\n')[0]
        pool_args.append({'sample_id_to_bam_file_path': sample_id_to_bam_file_path,
                          't_gene': t_gene,
                          't_gene_type': t_gene_type,
                          't_chrom': t_chrom,
                          't_start': int(t_start),
                          't_stop': int(t_stop),
                          'verbose': verbose,
                          'output_dir': output_dir})

    pool.map(find_splice_junctions_for_gene, pool_args)
    pool.close()
    pool.join()


def main():
    """If output directory is not provided, all files and intermediate files will be placed in a newly generated folder
    called sjd_output in the current working directory."""
    parser = argparse.ArgumentParser(description='Discover splice junctions from a list of bam files')
    parser.add_argument('transcript_file', metavar='transcript_file', type=str,
                        default='reference/gencode.comprehensive.splice.junctions.txt')
    parser.add_argument('--bam_folder', metavar='bam_folder', type=str)
    parser.add_argument('--bam_file', metavar='bam_file', type=str)
    parser.add_argument('--sample_id', metavar='sample_id', type=str, default='#')
    parser.add_argument('--threads', metavar='threads', type=int, default=10)
    parser.add_argument('--output_dir', metavar='output_dir', type=str)
    parser.add_argument('-v', action='store_true')
    parser.add_argument('-keep_gene_files', action='store_true')

    args = parser.parse_args()

    transcript_file = args.transcript_file
    bam_folder = args.bam_folder
    output_dir = args.output_dir
    threads = args.threads
    verbose = args.v

    # This is for the case where we want to process splice junctions one file at a time
    bam_file = args.bam_file
    sample_id = args.sample_id

    if not output_dir:
        output_dir = 'sjd_output'
        subprocess.call(['mkdir', output_dir])

    final_filename = '{}/final.txt'.format(output_dir)
    if os.path.isfile(final_filename):
        sys.exit('File already exists: {}. Please delete to rerun.'.format(final_filename))

    # For smaller sample sets, we can run them all at the same time
    if bam_folder:
        bam_file_paths = get_bam_files_in_folder(bam_folder)
        sample_id_to_bam_file_path = sample_ids = {
            get_id_from_bam_name(bam_file_path): bam_file_path for bam_file_path in bam_file_paths
        }

    # For larger sample sets, it is too time-consuming to run them all at the same time, so we do it one at a time and
    # then combine the results downstream
    elif bam_file:
        sample_id_to_bam_file_path = sample_ids = {
            sample_id: bam_file
        }

    # Mapping -- find the splice junctions across all samples, one gene per thread. Each thread generates a file.
    sys.stdout.write('>> Discovering splice junctions across all genes\n')
    map_splice_junction_discovery_across_genes(threads, transcript_file, sample_id_to_bam_file_path, verbose, output_dir)

    # Reducing -- combine all of the files generated into one megafile which includes all lines.
    sys.stdout.write('>> Combining all data into final file: {}\n'.format(final_filename))
    for root, dirs, files in os.walk(output_dir):
        open(final_filename, "w")
        subprocess.call("echo 'Gene\tType\tChrom\tStart\tEnd\tNTimesSeen\tNSamplesSeen\t{}' >> {}".format(
            '\t'.join(sample_ids),
            final_filename), shell=True)
        for file_counter, file in enumerate(sorted([f for f in files if f != final_filename])):
            gene_output_file = '{}/{}'.format(root, file)
            sys.stdout.write('\t* ({}/{}) adding data from {}\n'.format(file_counter+1, len(files), gene_output_file))
            subprocess.call('cat {} >> {}'.format(gene_output_file, final_filename), shell=True)
            # After feeding the contents into the aggregate file, delete the gene-level file
            if not args.keep_gene_files:
                subprocess.call('rm {}'.format(gene_output_file), shell=True)
        sys.stdout.write(">> Done! Output is in final file: {}\n".format(final_filename))


if __name__ == '__main__':
    main()


"""
Example run using folder of bams:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
--threads
2
--output_dir
/Users/erickofman/Documents/Projects/SpliceJunctions/PythonRNASeqJunctions/sample_id_adding
/Users/erickofman/Documents/Projects/SpliceJunctions/PythonRNASeqJunctions/reference/example.gene.NEB.list
--bam_folder
/Users/erickofman/Documents/Projects/SpliceJunctions/PythonRNASeqJunctions/example_bams


Example run using individual bam:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
--threads
2
--output_dir
/Users/erickofman/Documents/Projects/SpliceJunctions/PythonRNASeqJunctions/individual_bam
/Users/erickofman/Documents/Projects/SpliceJunctions/PythonRNASeqJunctions/reference/example.gene.NEB.list
--bam_file
/Users/erickofman/Documents/Projects/SpliceJunctions/PythonRNASeqJunctions/example_bams/Patient.D1.small.sorted.deduped.bam
--sample_id
D1
"""