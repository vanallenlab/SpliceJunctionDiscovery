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
    pre_N = cigar_string.split('M')[0]
    return [int(regex_num) for regex_num in re.findall('\d+', pre_N)][-1]


def write_cigar_debugging_info(cigar_string, num_matches_before_start, intron_length, pos, intron_start, intron_end):
    sys.stdout.write('>> {} ({} match, {} intronic)\t'
                     'read_start: {}\tintron_start: {}\tintron_end: {}\n'
                     .format(cigar_string,
                             num_matches_before_start,
                             intron_length,
                             pos,
                             intron_start,
                             intron_end))


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
                        intron_length = length_of_first_intronic_section(cigar_string)
                        num_matches_before_intron_start = num_of_matches_before_first_intronic_section(cigar_string)
                        intron_start_position = pos + num_matches_before_intron_start - 1
                        intron_end_position = pos + num_matches_before_intron_start + intron_length

                        if verbose:
                            # For debugging purposes, if verbose is true, output some intermediate information
                            write_cigar_debugging_info(cigar_string, num_matches_before_intron_start, intron_length,
                                                       pos, intron_start_position, intron_end_position)

                        unique_splice_junction = '{},{},{}'.format(t_chrom, intron_start_position, intron_end_position)
                        splice_junctions[unique_splice_junction] += 1

                    # Account for very small retained introns (they must be small for the whole intron to have been
                    # captured within the length of one read), e.g. 13M221400N34M2658N29M
                    if cigar_string.count('N') == 2:
                        section_after_first_intron = cigar_string.split('N')[1]
                        num_matches_before_second_intron = \
                            num_of_matches_before_first_intronic_section(section_after_first_intron)
                        second_intron_length = length_of_first_intronic_section(section_after_first_intron)
                        second_intron_start_position = intron_end_position + num_matches_before_second_intron
                        second_intron_end_position = second_intron_start_position + second_intron_length

                        # For debugging purposes, if verbose is true, output some intermediate information
                        if verbose:
                            write_cigar_debugging_info(cigar_string, num_matches_before_second_intron,
                                                       second_intron_length, pos, second_intron_start_position,
                                                       second_intron_end_position)

                        unique_splice_junction = '{},{},{}'.format(t_chrom, second_intron_start_position,
                                                                   second_intron_start_position)
                        splice_junctions[unique_splice_junction] += 1

        except IndexError:
            pass

    return splice_junctions


def get_id_from_bam_name(bam_file_path):
    """Strip away the .bam suffix and return the sample ID"""
    filename, file_extension = os.path.splitext(os.path.basename(bam_file_path))
    return filename


def summarize_splice_junctions(gene, global_event_counts, output_dir='.'):
    """Given a dictionary summarizing the global and per sample occurrence counts of unique splice junctions,
    output a file summarizing the results"""
    output_file_name = '{}/{}_{}_junctions.txt'.format(output_dir, gene, len(global_event_counts))
    with open(output_file_name, 'w') as f:
        ordered_events = sorted(global_event_counts.keys())
        for event in ordered_events:
            per_sample_counts = global_event_counts[event]
            gene, gene_type, chrom, start, end = event.split(',')

            total_count = 0
            sample_summaries = []
            for sample, count in per_sample_counts.items():
                total_count += count
                sample_summaries.append('{}:{}'.format(sample, count))

            sample_summary_string = '{}'.format(','.join(sample_summaries))
            num_samples_with_this_event = len(per_sample_counts)

            entry = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(gene, gene_type, chrom, start, end, total_count,
                                                               num_samples_with_this_event, sample_summary_string)
            f.write(entry)


def find_splice_junctions_for_gene(pool_arguments):
    """Find the splice junctions across all samples for one gene, writing the results to a file
    specific to this gene"""
    pa = pool_arguments
    # A dictionary to keep track of, for each unique splice junction event, how many times it appears in each of the
    # samples being considered as part of this analysis for this gene.
    global_event_counts = defaultdict(lambda: defaultdict())
    for bam_file_path, sample_id in zip(pa.get('bam_file_paths'), pa.get('sample_ids')):
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
    summarize_splice_junctions(pa.get('t_gene'), global_event_counts, pa.get('output_dir'))


def map_splice_junction_discovery_across_genes(threads, transcript_file, bam_paths, sample_ids, verbose=False, output_dir=None):
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
        pool_args.append({'bam_file_paths': bam_paths,
                          'sample_ids': sample_ids,
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
    parser.add_argument('bam_folder', metavar='bam_folder', type=str, default='bams')
    parser.add_argument('--threads', metavar='threads', type=int, default=10)
    parser.add_argument('--output_dir', metavar='output_dir', type=str)
    parser.add_argument('-v', action='store_true')
    parser.add_argument('-keep_gene_files', action='store_true')

    args = parser.parse_args()

    transcript_file = args.transcript_file
    bam_folder = args.bam_folder
    output_dir = args.output_dir
    threads = args.threads

    if not output_dir:
        output_dir = 'sjd_output'
        subprocess.call(['mkdir', output_dir])

    final_filename = '{}/final.txt'.format(output_dir)
    if os.path.isfile(final_filename):
        sys.exit('File already exists: {}. Please delete to rerun.'.format(final_filename))

    verbose = args.v
    bam_file_paths = get_bam_files_in_folder(bam_folder)
    sample_ids = [get_id_from_bam_name(bam_file_path) for bam_file_path in bam_file_paths]

    # Mapping -- find the splice junctions across all samples, one gene per thread. Each thread generates a file.
    sys.stdout.write('>> Discovering splice junctions across all genes\n')
    map_splice_junction_discovery_across_genes(threads, transcript_file, bam_file_paths, sample_ids, verbose, output_dir)

    # Reducing -- combine all of the files generated into one megafile which includes all lines.
    sys.stdout.write('>> Combining all data into final file: {}\n'.format(final_filename))
    for root, dirs, files in os.walk(output_dir):
        open(final_filename, "w")
        subprocess.call("echo 'Gene\tType\tChrom\tStart\tEnd\tNTimesSeen\tNSamplesSeen\tSamples:NSeen' >> {}".format(
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
