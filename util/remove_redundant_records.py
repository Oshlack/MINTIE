import argparse
import pandas as pd
import re
import sys
import fileinput
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument(dest='fasta_file')
args = parser.parse_args()

fasta_file = args.fasta_file
with fileinput.input() as handle:
    contigs = []
    for record in SeqIO.parse(handle, 'fasta'):
        record.description = record.id
        if record.id not in contigs:
            sys.stdout.write(record.format('fasta'))
            contigs.append(record.id)
