import argparse
import pandas as pd
import sys


def main(args):
    normalized_file = args.normalized
    df = pd.read_csv(normalized_file, sep='\t', header=0)

    sys.stdout.write("{} total splice junctions\n".format(len(df)))
    sys.stdout.write("Filtering out {} annotated splice junctions...\n".format(len(df[df['SeenInTranscriptsList'] == 1])))
    df = df[df['SeenInTranscriptsList'] == 0]
    sys.stdout.write("{} non-canonical splice junctions\n".format(len(df)))
    print("hello")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''Filter normalized splice junctions file''')
    parser.add_argument('-normalized', help='normalized splice junctions file', type=str)
    args = parser.parse_args()
    main(args)
