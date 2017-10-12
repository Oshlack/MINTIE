import pandas as pd
import argparse
import os
import sys
import numpy as np
from intervaltree import Interval, IntervalTree

parser = argparse.ArgumentParser()
parser.add_argument(dest='annotation')
parser.add_argument(dest='assembly')
parser.add_argument(dest='partition')

args = parser.parse_args()
part_path = args.partition
ann_path = args.annotation
ass_path = args.assembly

bed_paths = [part_path, ann_path, ass_path]
if not all([os.path.exists(path) for path in bed_paths]):
    sys.exit('Not all bed files supplied exist!')

partitioned = pd.read_csv(part_path, sep='\t', header=None)
annotation = pd.read_csv(ann_path, sep='\t', header=None)
assembled = pd.read_csv(ass_path, sep='\t', header=None)

def make_tree(bed_ranges):
    tree = IntervalTree()
    for idx, row in bed_ranges.iterrows():
        start, end = row[1], row[2]
        tree[start:end] = '%d-%d' % (start, end)
    return(tree)

contigs = np.unique(partitioned[0].values)

for contig in contigs:
    part_contig = partitioned[partitioned[0].values==contig]
    ann_contig = annotation[annotation[0].values==contig]
    ass_contig = assembled[assembled[0].values==contig]
    
    part_tree = make_tree(part_contig)
    ann_tree = make_tree(ann_contig)
    ass_tree = make_tree(ass_contig)

    presence = []
    for interval in part_tree:
        outstr = ''
        if ann_tree.overlaps(interval):
            outstr = 'Annotation'
        if ass_tree.overlaps(interval):
            outstr = 'Assembly' if outstr == '' else '%s,%s' % (outstr, 'Assembly')
        presence.append(outstr)
    
    part_contig = part_contig.assign(presence = pd.Series(presence, index=part_contig.index))
    part_contig.to_csv(sys.stdout, header=None, columns=None, index=False, sep='\t')

