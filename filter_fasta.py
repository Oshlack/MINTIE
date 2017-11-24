########################################################
# Author: Marek Cmero
# Take list of transcripts/contigs and filter out any
# sequences not in this list from target fasta file
########################################################

import argparse
import pandas as pd
import re
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument(dest='fasta_file')
parser.add_argument(dest='tx_list')
args = parser.parse_args()

tx_list_file = args.tx_list
fasta_file = args.fasta_file

tx_list = pd.read_csv(tx_list_file, sep='\t', header=None)

handle = open(fasta_file, 'r')
for record in SeqIO.parse(handle, 'fasta'):
    if record.id in tx_list[0].values:
        print(record.format('fasta'))
handle.close()
