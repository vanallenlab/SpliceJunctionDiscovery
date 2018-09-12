import argparse
import pandas as pd
import sys
from os import path


def main(args):
    normalized_file = args.normalized
    min_supporting_reads = args.min_supporting_reads
    output_dir = args.output_dir
    df = pd.read_csv(normalized_file, sep='\t', header=0)

    sys.stdout.write("{} total splice junctions\n".format(len(df)))
    sys.stdout.write("Filtering out {} annotated splice junctions...\n".format(len(df[df['SeenInTranscriptsList'] == 1])))
    df = df[df['SeenInTranscriptsList'] == 0]
    sys.stdout.write("{} non-canonical splice junctions\n".format(len(df)))

    sys.stdout.write("Filtering out sites where all samples had below {} reads...\n".format(min_supporting_reads))

    start_position_of_sample_count_columsn = 7
    end_position_of_sample_count_columns = list(df.columns).index('Annotations')
    sample_names = df.columns[start_position_of_sample_count_columsn: end_position_of_sample_count_columns]

    eval_command = '&'.join(
        ["((df[\'{}\'] >= {}) | (df[\'{}\'] == 0))".format(sample_name, min_supporting_reads, sample_name)
         for sample_name in sample_names])
    print('Command: {}\n'.format(eval_command))

    and_sub_clauses_for_df = eval(eval_command)

    candidate_sites = df[and_sub_clauses_for_df]
    sys.stdout.write('{} candidate sites\n'.format(len(candidate_sites)))

    candidate_sites.to_csv('{}/candidates.{}'.format(output_dir, path.basename(normalized_file)), sep='\t')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''Filter normalized splice junctions file''')
    parser.add_argument('-normalized', help='normalized splice junctions file', type=str)
    parser.add_argument('-min_supporting_reads', help='min supporting number of reads', type=int, default=5)
    parser.add_argument('-output_dir', help='Where to output the output file', type=str, default='./')
    args = parser.parse_args()
    main(args)
