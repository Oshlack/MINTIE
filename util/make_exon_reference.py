###########################################################
# Author: Marek Cmero
# Takes a GTF annotation file of transcripts (e.g. gencode)
# and outputs one line per transcript with comma separated
# lists of exon starts and ends
###########################################################

#TODO: clean-up/document this script better

import argparse
import pandas as pd
import numpy as np
import os
from pybedtools import BedTool

parser = argparse.ArgumentParser()
parser.add_argument(dest='tx_info')

args = parser.parse_args()
tx_info = args.tx_info

print('Reading in gene reference...')
genref = pd.read_csv(tx_info, sep='\t', comment='#', header=None, low_memory=False)
pd.options.mode.chained_assignment = None #to get rid of those pesky annoying pandas warnings

def featuretype_filter(feature, featuretype):
    '''
    from http://daler.github.io/pybedtools/3-brief-examples.html
    Only passes features with the specified *featuretype*
    '''
    if feature[2] == featuretype:
        return True
    return False

def subset_featuretypes(g, featuretype):
    '''
    from http://daler.github.io/pybedtools/3-brief-examples.html
    Returns the filename containing only `featuretype` features.
    '''
    return g.filter(featuretype_filter, featuretype).saveas().fn

def get_attribute(g, attribute):
    genes = []
    for feature in g:
        try:
            genes.append(feature[attribute])
        except AttributeError:
            genes.append('')
    return genes

print('loading bedfile and extracting exons...')
g = BedTool(tx_info)
exons = BedTool(subset_featuretypes(g, 'exon'))

print('validating and sorting records...')
exons = exons.remove_invalid().sort()

print('extracting attributes...')
exon_pd = pd.DataFrame([(e['chrom'], e['start'], e['end'], e['strand']) for e in exons],
                        columns=['chrom', 'exonStarts', 'exonEnds', 'strand'])
exon_pd['exonStarts'] = exon_pd['exonStarts'].map(int)
exon_pd['exonEnds'] = exon_pd['exonEnds'].map(int)
exon_pd['transcript'] = get_attribute(exons, 'transcript_id')
exon_pd['gene'] = get_attribute(exons, 'gene_name')
exon_pd = exon_pd[exon_pd.gene != '']

print('building exon reference...')
aggregator_sense = {'exonStarts': lambda x: ','.join(x.map(str)),
                    'exonEnds': lambda x: ','.join(x.map(str))}
all_exons = exon_pd.groupby(['transcript', 'chrom', 'gene'], as_index=False, sort=False).agg(aggregator_sense)

print('writing output...')
outdir = os.path.dirname(tx_info)
outdir = '.' if outdir == '' else outdir
outfile = '%s/%s.info' % (outdir, '.'.join(os.path.basename(tx_info).split('.')[:-1]))
all_exons = all_exons[['transcript', 'chrom', 'exonStarts', 'exonEnds', 'gene']]
all_exons.to_csv(outfile, sep='\t', index=False)
