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

ec_files = args.inputs[:-1]
outfile = args.inputs[-1]

def load_ecs(ec_file):
    '''
    Read Salmon equivalence class output file
    and return dictionary with all transcripts,
    equivalence class transcript IDs and counts
    '''
    ec_df = pd.DataFrame([line.strip().split('\t') for line in open(ec_file, 'r')], columns=None)
    ec_rows = np.array([bool(re.search('^\d+$', val)) for val in ec_df[0].values])
    ecs = ec_df[ec_rows][2:]

    # extract counts and transcript IDs
    transcripts = ec_df[np.invert(ec_rows)][0].values
    ec_vals = ecs.apply(lambda x: [int(val) for val in x if val], axis=1)
    counts = ec_vals.apply(lambda x: x[-1]).values
    tx_ids = ec_vals.apply(lambda x: x[1:-1]).values

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

cancer = load_ecs(ec_files[0])
controls = [load_ecs(file) for file in ec_files[1:]]

# build EC dictionary
ec_dict = build_ec_dict({}, cancer, 'cancer')
for idx, control in enumerate(controls):
    control_name = 'control%d' % (idx+1)
    ec_dict = build_ec_dict(ec_dict, control, control_name)

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
counts['transcript'] = np.array([cancer['transcripts'][tidx] for tidx in counts.tx_id.map(int).values])

counts.to_csv(outfile, sep='\t', index=False)
