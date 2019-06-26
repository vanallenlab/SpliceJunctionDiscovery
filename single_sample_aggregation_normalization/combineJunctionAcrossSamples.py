import argparse 
import os 
from itertools import chain
from glob import glob
import pandas as pd 
import sys
import json 
from collections import OrderedDict, defaultdict 
import numpy as np

def main(args):

    sample_ids = []

    all_juncs = []
    all_juncs_dict = OrderedDict()

    for f in glob('{}/*.splicejunctions.reformatted.tsv'.format(args.input_junction_folder)):
        splice_df = pd.read_csv(f, sep='\t')

        samptimes_list = [json.loads(k) for k in splice_df.samptimes.values]
        sample_id = np.unique([k for d in samptimes_list for k in d.keys()])[0]
        
        sys.stdout.write('Sample ID {}\n'.format(sample_id))

        sample_ids.append(sample_id)
        
    for f in glob('{}/*.splicejunctions.reformatted.tsv'.format(args.input_junction_folder)):

        splice_df = pd.read_csv(f, sep='\t')
        splice_df['samptimes'] = splice_df['samptimes'].apply(lambda x: json.loads(x))
        splice_junctions = splice_df.to_dict('records')

        for junction in splice_junctions:
            chrom, start, stop, samptimes = junction.get('chrom'), junction.get('start'), junction.get('stop'), junction.get('samptimes')

            (samp_id, ntimes), = samptimes.items()
            junction_string = '{}:{}-{}'.format(chrom, start, stop)

 
            if junction_string in all_juncs_dict:
                if samp_id in all_juncs_dict[junction_string]['samptimes']:
                    all_juncs_dict[junction_string]['samptimes'][samp_id] += ntimes
                else:
                    all_juncs_dict[junction_string]['samptimes'][samp_id] = ntimes
                all_juncs_dict[junction_string]['ntimes'] += ntimes
                all_juncs_dict[junction_string]['nsamp'] += 1
                

            else:
                all_juncs_dict[junction_string] = ({
                    'gene': junction.get('gene'),
                    'gene_type': junction.get('gene_type'),
                    'chrom': chrom,
                    'start': start,
                    'stop': stop,
                    'ntimes': ntimes,
                    'nsamp': 1,
                })

                all_juncs_dict[junction_string]['samptimes'] = defaultdict(int)
                all_juncs_dict[junction_string]['samptimes'][samp_id] = ntimes

        
        
        #data = pd.read_csv(f, sep='\t')
        #all_juncs.append(data)

    

    #df_all = pd.concat(all_juncs)

    junc_values = list(all_juncs_dict.values())
    new_junc_values = []
    for juncs in junc_values:
        temp_dict = juncs
        temp_dict['nsamp'] = np.count_nonzero(np.asarray(list(temp_dict['samptimes'].values())))
        temp_dict['samptimes'] = json.dumps(temp_dict['samptimes'])
        new_junc_values.append(temp_dict)


    output_df = pd.DataFrame(new_junc_values)
    output_file_name = '{}/{}_allsplicejunctions.tsv'.format(args.output_folder, args.sample_set_id)

    output_df.to_csv(output_file_name, sep='\t', index=False)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-input_junction_folder', help='folder containing all splice junctions from all samples')
    parser.add_argument('-output_folder', type=str, default='outputs', help='folder for outputs')
    parser.add_argument('-sample_set_id', type=str, default='sample_set', help='sample_set_id')
    args = parser.parse_args()

    main(args)