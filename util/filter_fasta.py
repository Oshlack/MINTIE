########################################################
# Author: Marek Cmero
# Take list of transcripts/contigs and filter out any
# sequences not in this list from target fasta file
########################################################

import argparse
import pandas as pd
import re
import sys
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument(dest='fasta_file')
parser.add_argument(dest='tx_list')
parser.add_argument('--col_id', dest='col_id', default=0)
args = parser.parse_args()

tx_list_file = args.tx_list
fasta_file = args.fasta_file
col_id = args.col_id

header = None
try:
    # try to get col numeric index, if this fails, assume the file is headered
    col_id = int(col_id)
except ValueError:
    header = 0

tx_list = pd.read_csv(tx_list_file, sep='\t', header=header)

lookup_list = set(tx_list[col_id].values.tolist())

handle = open(fasta_file, 'r')
for record in SeqIO.parse(handle, 'fasta'):
    record.description = record.id
    if record.id in lookup_list:
        sys.stdout.write(record.format('fasta'))
handle.close()
