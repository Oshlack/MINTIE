###########################################################
# Author: Marek Cmero
# Takes a GTF annotation file of transcripts (e.g. gencode)
# and outputs one line per transcript with comma separated
# lists of exon starts and ends
###########################################################

import argparse
import pandas as pd
import numpy as np
import os

parser = argparse.ArgumentParser()
parser.add_argument(dest='tx_info')

args = parser.parse_args()
tx_info = args.tx_info

genref = pd.read_csv(tx_info, sep='\t', comment='#', header=None)
pd.options.mode.chained_assignment = None #to get rid of those pesky annoying pandas warnings

exons = genref[genref[2] == 'exon']
exons[3] = exons[3]-1 # minus one for each start coordinate to make them 0-based
exons[9] = [info.split(';')[1].strip().split(' ')[1].replace('"','') for info in genref[genref[2] == 'exon'][8]]
exons[10] = [info.split(';')[4].strip().split(' ')[1].replace('"','') for info in genref[genref[2] == 'exon'][8]]

# join starts and ends together
exons_sense = exons[exons[6] == '+']
aggregator_sense = {3: lambda x: ','.join(x.map(str)),
                    4: lambda x: ','.join(x.map(str))}
exons_sense = exons_sense.groupby([9,0,10], as_index=False, sort=False).agg(aggregator_sense)

# same thing but reverse the exon order
exons_antisense = exons[exons[6] == '-']
aggregator_antisense = {3: lambda x: ','.join(x[::-1].map(str)),
                        4: lambda x: ','.join(x[::-1].map(str))}
exons_antisense = exons_antisense.groupby([9,0,10], as_index=False, sort=False).agg(aggregator_antisense)

all_exons = exons_sense.append(exons_antisense, ignore_index=True)
all_exons = all_exons[[9, 0, 3, 4, 10]]
all_exons.columns = ['transcript', 'chrom', 'exonStarts', 'exonEnds', 'gene']

outfile = '%s/%s.info' % (os.path.dirname(tx_info), '.'.join(os.path.basename(tx_info).split('.')[:-1]))
all_exons.to_csv(outfile, sep='\t', index=False)
