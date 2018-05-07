import sys
import argparse
import pandas as pd
import glob
import os
import subprocess
import re
from collections import defaultdict

SAM_CHROM_COL_INDEX = 2
SAM_CHROM_POS_INDEX = 3
SAM_CIGAR_COL_INDEX = 5
SAM_COL_INDICES = [SAM_CHROM_COL_INDEX, SAM_CHROM_POS_INDEX, SAM_CIGAR_COL_INDEX]


def get_group_name(transcript_file):
    return 'Group'


def get_bam_files_in_folder(bam_folder):
    """Get a list of all the bam file paths in a given folder"""
    bam_file_paths = []
    for root, dirs, files in os.walk(bam_folder):
        for file in files:
            if file.endswith(".bam"):
                 bam_file_paths.append(os.path.join(root, file))
    return bam_file_paths


def length_of_first_intronic_section(cigar_string):
    """The intron length will be the last regex integer to match the string before the 'N' character"""
    return int(re.findall('\d+', cigar_string.split('N')[0])[-1])


def num_of_matches_before_first_intronic_section(cigar_string):
    pre_N = cigar_string.split('M')[0]
    return [int(regex_num) for regex_num in re.findall('\d+', pre_N)][-1]


def find_splice_junctions(bam_file_path, t_chrom, t_start, t_stop, verbose=False):
    t_location = '{}:{}-{}'.format(t_chrom, t_start, t_stop)
    # Samtools view options:
    # ~~~~~~~~~~~~~~~~~~~~~~
    # -F INT Do not output alignments with any bits set in INT present in the FLAG field
    # -q INT    Skip alignments with MAPQ smaller than INT
    pos = 'chr{}:{}-{}'.format(t_chrom, t_start, t_stop)
    samtools_view_command = 'samtools view -F 256 -F 512 -F 1024 -F 2048 -q 60 {} {} {}'.format(bam_file_path,
                                                                                                pos,
                                                                                                t_location)
    # Run the samtools view command using subprocess, and place the stdout output into the sam_view variable as a string
    sam_view = subprocess.check_output(samtools_view_command.split())

    # Split the sam_view output string into each line comprising it
    sam_lines = str(sam_view).split('\\n')

    # A dictionary that will keep track of how many times a particular junction is found
    splice_junctions = defaultdict(int)

    # Only look at the exon junction spanning lines, for which the CIGAR string will contain an 'N' character
    for line in sam_lines:
        try:
            split_line = line.split('\\t')
            cigar_string = split_line[SAM_CIGAR_COL_INDEX]
            pos = int(split_line[SAM_CHROM_POS_INDEX])
            if 'N' in cigar_string and t_start <= pos <= t_stop:
                    if cigar_string.count('N') <= 2:
                        # Account for only one intronic region, or only the first of two in the case there are two
                        intron_length = length_of_first_intronic_section(cigar_string)
                        num_matches_before_intron_start = num_of_matches_before_first_intronic_section(cigar_string)
                        intron_start_position = pos + num_matches_before_intron_start - 1
                        intron_end_position = pos + num_matches_before_intron_start + intron_length

                        if verbose:
                            # For debugging purposes, if verbose is true, output some intermediate information
                            sys.stdout.write('{} ({} match, {} intronic)\t'
                                             'read_start: {}\tintron_start: {}\tintron_end: {}\n'
                                             .format(cigar_string,
                                                     num_matches_before_intron_start,
                                                     intron_length,
                                                     pos,
                                                     intron_start_position,
                                                     intron_end_position))

                        unique_splice_junction = '{},{},{}'.format(t_chrom, intron_start_position, intron_end_position)
                        splice_junctions[unique_splice_junction] += 1

                    elif cigar_string.count('N') == 2:
                        # Account for very small retained introns (they must be small for the whole intron to have been
                        # captured within the length of one read)
                        pass
        except IndexError:
            pass

    return splice_junctions


def get_id_from_bam_name(bam_file_path):
    """Strip away the .bam suffix and return the sample ID"""
    filename, file_extension = os.path.splitext(os.path.basename(bam_file_path))
    return filename


def summarize_splice_junctions(global_event_counts, output_file_name='PythonSpliceJunctionSummary.txt'):
    """Given a dictionary summarizing the global and per sample occurrence counts of unique splice junctions,
    output a file summarizing the results"""
    with open(output_file_name, 'w') as f:
        f.write('gene\ttype\tchrom\tstart\tend\tinstances\tnum_samples\tsample_level_info\n')
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


def main():
    parser = argparse.ArgumentParser(description='Find splice junction sites')
    parser.add_argument('transcript_file', metavar='transcript_file', type=str)
    parser.add_argument('bam_folder', metavar='bam_folder', type=str)
    parser.add_argument('--output_dir', metavar='output_dir', type=str)
    parser.add_argument('-v', action='store_true')

    args = parser.parse_args()

    transcript_file = args.transcript_file
    bam_folder = args.bam_folder
    output_dir = args.output_dir
    verbose = args.v

    group_name = get_group_name(transcript_file)

    bam_file_paths = get_bam_files_in_folder(bam_folder)
    print(bam_file_paths)

    # A dictionary to keep track of, for each unique splice junction event, how many times it appears in each of the
    # samples being considered as part of this analysis.
    global_event_counts = defaultdict(lambda: defaultdict())

    transcript_lines = open(transcript_file).readlines()
    for t in transcript_lines:
        t_gene, t_enst, _, t_chrom, t_start, t_stop, t_gene_type = t.split('\t')
        t_gene = t_gene.split('\n')[0]
        t_gene_type = t_gene_type.split('\n')[0]
        for bam_file_path in bam_file_paths:
            patient_id = get_id_from_bam_name(bam_file_path)
            splice_junctions = find_splice_junctions(bam_file_path,
                                                     int(t_chrom),
                                                     int(t_start),
                                                     int(t_stop),
                                                     verbose=verbose)
            for splice_junction, count in splice_junctions.items():
                splice_junction_with_gene_info = '{},{},{}'.format(t_gene, t_gene_type, splice_junction)
                global_event_counts[splice_junction_with_gene_info][patient_id] = count

    summarize_splice_junctions(global_event_counts)


if __name__ == '__main__':
    main()