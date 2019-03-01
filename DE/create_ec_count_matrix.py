########################################################
# Author: Marek Cmero
# Match salmon equivalence class output files from multiple
# samples and create a matched transcript and EC matrix
########################################################

import os
import argparse
import re
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('inputs', nargs='*')
args = parser.parse_args()

ec_files = args.inputs[:-2]
sample_names = args.inputs[-2]
outfile = args.inputs[-1]

def load_ecs(ec_file):
    '''
    Read Salmon equivalence class output file
    and return dictionary with all transcripts,
    equivalence class transcript IDs and counts
    '''
    ec_df = pd.read_csv(ec_file, header=None)
    ec_df = ec_df[0].apply(lambda x: x.split('\t'))
    transcripts = ec_df[ec_df.apply(len)==1][2:]
    transcripts = [t for tx in transcripts for t in tx]

    # extract counts and transcript IDs
    ec_df = ec_df[ec_df.apply(len)>1]
    counts = ec_df.apply(lambda x: int(x[-1])).values
    tx_ids = ec_df.apply(lambda x: x[1:-1]).values

    output = {'transcripts': transcripts,
              'tx_ids': tx_ids,
              'counts': counts}

    return(output)

def build_ec_dict(ec_dict, sample, name):
    '''
    Build equivalence class dictionary of
    counts using transcript IDs as keys
    '''
    tx_ids = map(lambda x: '|'.join(list(map(str, x))), sample['tx_ids'])
    for tx_id, count in zip(tx_ids, sample['counts']):
        if tx_id not in ec_dict.keys():
            ec_dict[tx_id] = {}
        ec_dict[tx_id][name] = count
    return(ec_dict)

sample_ecs = [load_ecs(file) for file in ec_files]
sample_names = sample_names.split(',')

# build EC dictionary
ec_dict = {}
for idx, sample_ec in enumerate(sample_ecs):
    ec_dict = build_ec_dict(ec_dict, sample_ec, sample_names[idx])

# construct counts dataframe
counts = pd.DataFrame(ec_dict).transpose().fillna(0)
ec_names = ['ec%d' % (i+1) for i in range(len(counts))]
counts = counts.assign(ec_names=ec_names)
counts = counts.assign(tx_ids=counts.index.values)

# split ECs with multiple transcript IDs into separate rows
tmp = counts.tx_ids.str.split('|').apply(pd.Series, 1).stack()
tmp.index = tmp.index.droplevel(-1)
tmp = pd.DataFrame(tmp, columns=['tx_id'])
tmp.name = 'tx_ids'
counts = counts.join(tmp)
del counts['tx_ids']

# transcript IDs > names
counts['transcript'] = np.array([sample_ecs[0]['transcripts'][tidx] for tidx in counts.tx_id.map(int).values])

counts.to_csv(outfile, sep='\t', index=False)
