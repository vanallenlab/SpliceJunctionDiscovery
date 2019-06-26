import argparse 
import os 
import pandas as pd
import json 
import sys

def junctions_from_file(splice_file):

    junctions = []

    with open(splice_file, 'r') as f:
        sample_id = [k for k in f.readline().split("\t")[7:]][0].rstrip()
        file_lines = f.readlines()[1:]

        for l in file_lines:
            gene, gene_type, chrom, start, stop, ntimes, nsamp = l.strip().split("\t")[0:7]
            sample_counts = [int(c.strip()) for c in l.split("\t")[7:]]

            junctions.append({
                'gene': gene,
                'gene_type': gene_type,
                'chrom': chrom,
                'start': start,
                'stop': stop,
                'ntimes': ntimes,
                'nsamp': nsamp,
                'samptimes': json.dumps({sample_id: v for v in sample_counts})
            })

    return junctions


def main(args):

    junc = junctions_from_file(args.junction_file)
    df = pd.DataFrame(junc)

    junction_file_base = os.path.basename(args.junction_file)

    output_file_name = '{}/{}.splicejunctions.reformatted.tsv'.format(args.output_folder, junction_file_base)

    df.to_csv(output_file_name, sep='\t', index=False)





if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-junction_file', help='.txt file containing output from SpliceJunctionDiscovery')
    parser.add_argument('-output_folder', type=str, default='outputs', help='folder for outputs')
    args = parser.parse_args()

    main(args)