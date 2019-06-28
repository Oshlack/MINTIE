'''
Module      : create_ec_count_matrix
Description : Perform contig annotation.
Copyright   : (c) Marek Cmero, Dec 2018
License     : MIT
Maintainer  : MAREK.CMERO@MCRI.EDU.AU
Portability : POSIX
Match salmon equivalence class output files from multiple
samples and create a matched transcript and EC matrix
'''
import pandas as pd
import numpy as np
from argparse import ArgumentParser

def parse_args():
    '''
    Parse command line arguments.
    Returns Options object with command line argument values as attributes.
    Will exit the program on a command line error.
    '''
    description = 'Create EC counts matrix from Salmon input'
    parser = ArgumentParser(description=description)
    parser.add_argument('inputs', nargs='*')

    return parser.parse_args()


def load_ecs(ec_file):
    '''
    Read Salmon equivalence class output file
    and return dictionary with all transcripts,
    equivalence class transcript IDs and counts
    '''
    ec_df = pd.read_csv(ec_file, header=None)
    ecc = [x.split('\t') for x in ec_df[0].values if len(x.split('\t')) > 1]
    ec_txs = ['|'.join(x[1:-1]) for x in ecc]
    ec_counts = [int(x[-1]) for x in ecc]

    return ec_txs, ec_counts

def build_ec_matrix(sample_ecs, sample_names):
    '''
    Build equivalence class dictionary of
    counts using transcript IDs as keys
    '''
    ec_dicts = {}
    for idx,sec in enumerate(sample_ecs):
        ec_dict = {}
        ec_list = [ec for ec in zip(sec[0], sec[1])]
        ec_dict.update(ec_list)
        ec_dicts[sample_names[idx]] = ec_dict

    return pd.DataFrame(ec_dicts).fillna(0)

def construct_dataframe(ec_matrix, tx_lookup):
    '''
    - split transcript IDs
    - add EC names
    - convert transcript IDs to transcript names
    '''
    tmp = [(eid, etxs.split('|')) for eid, etxs in enumerate(ec_matrix.index.values)]
    tmp = [(tid, eid) for eid, tids in tmp for tid in tids]
    tx_ids, ec_ids = zip(*tmp)

    ec_names = ['ec%d' % (i+1) for i in ec_ids]
    ec_matrix = ec_matrix.iloc[np.array(ec_ids)]
    ec_matrix = ec_matrix.assign(ec_names=ec_names)
    ec_matrix = ec_matrix.assign(tx_ids=tx_ids)

    tx_names = [tx_lookup[tx_id] for tx_id in ec_matrix.tx_ids.map(int).values]
    ec_matrix = ec_matrix.assign(transcript=tx_names)

    return ec_matrix

def get_tx_lookup(ec_file):
    ec_df = pd.read_csv(ec_file, header=None)
    tx_lookup = [x for x in ec_df[0].values if len(x.split('\t')) == 1][2:]
    return np.array(tx_lookup)

def main():
    args = parse_args()
    ec_files = args.inputs[:-2]
    sample_names = args.inputs[-2]
    outfile = args.inputs[-1]

    print('Loading ECs...')
    sample_ecs = [load_ecs(file) for file in ec_files]
    sample_names = sample_names.split(',')

    print('Building EC matrix...')
    ec_matrix = build_ec_matrix(sample_ecs, sample_names)

    print('Constructing dataframe...')
    tx_lookup = get_tx_lookup(ec_files[0])
    ec_matrix = construct_dataframe(ec_matrix, tx_lookup)

    print('Writing output...')
    ec_matrix.to_csv(outfile, sep='\t', index=False)

if __name__ == '__main__':
    main()
