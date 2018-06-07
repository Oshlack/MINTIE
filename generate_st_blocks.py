#########################################################
# Author: Marek Cmero
#########################################################

import argparse
import pandas as pd
import numpy as np
import gffutils
from intervaltree import Interval, IntervalTree

parser = argparse.ArgumentParser()
parser.add_argument(dest='gff_file',
                    help='''GTF reference annotation''')
parser.add_argument(dest='out_bed',
                    help='''Output bed file of exonic coordinates.''')

args        = parser.parse_args()
gff_file    = args.gff_file
out_bed     = args.out_bed

db_name = '%s.db' % gff_file.split('.')[:-1][0]
db = gffutils.create_db(gff_file, db_name, force=True, disable_infer_transcripts=True, disable_infer_genes=True, verbose=True)

exons = {}
for gene in db.features_of_type('gene'):
    gene_name = gene.attributes['gene_name'][0]
    gene_intervals = IntervalTree()
    for exon in db.children(gene, featuretype='exon'):
        if exon.end - exon.start > 0:
            gene_intervals[exon.start:exon.end] = [exon.strand, exon.attributes['exon_number'][0]]

    IntervalTree.split_overlaps(gene_intervals)
    if gene.chrom not in exons.keys():
        exons[gene.chrom] = {}
    exons[gene.chrom][gene_name] = gene_intervals

exon_df = []
for chrom in exons:
    for gene in exons[chrom]:
        intervals = [[gene, chrom, x[0], x[1], x[2][1], '.', x[2][0]] for x in exons[chrom][gene]]
        exon_df.extend(intervals)

exon_df = pd.DataFrame(exon_df).drop_duplicates()
exon_df = exon_df.groupby([0,1,2,3,5,6])[4].apply(lambda x: '|'.join(x)).to_frame().reset_index()
exon_df.columns = ['gene', 'chrom', 'start', 'end', 'score', 'strand', 'exons']

exon_df['name'] = exon_df.gene + ' ' + exon_df.exons
exon_df = exon_df[['chrom', 'start', 'end', 'name', 'score', 'strand']]
exon_df.chrom = exon_df.chrom.map(lambda x: x.split('chr')[1])
exon_df.loc[exon_df.chrom == 'M', 'chrom'] = 'MT'

starts = exon_df.chrom.map(str) + ':' + exon_df.start.map(str)
ends = exon_df.chrom.map(str) + ':' + exon_df.end.map(str)
non_adjacent_blocks = np.array([s not in ends.values for s in starts])
exon_df.loc[non_adjacent_blocks, 'start'] = exon_df.start - 1 # only offset blocks that are not adjacent
exon_df.to_csv(out_bed, sep='\t', index=False, header=False)
