#########################################################
# Author: Marek Cmero
# Take an samfile of aligned contigs, and filter out any
# contigs that match annotated transcript characteristics
# i.e. - genomic gaps <7bp
#      - no soft or hard clips >30bp
#      - all junctions are known (in annotation)
#########################################################

import argparse
import pandas as pd
import numpy as np
import pysam
import os
import re
from Bio import SeqIO

# cutoff parameters
gap_min = 7
clip_min = 30

# CIGAR specification codes
gaps = {'deletion': 1, 'insertion': 2, 'silent_deletion': 6}
clips = {'soft': 4, 'hard': 5}

def get_juncs(tx):
    '''
    return list of junctions in form
    [(chr, start, end)] from transcript
    info file
    '''
    starts = tx['exonStarts'].split(',')[1:-1]
    ends = tx['exonEnds'].split(',')[:-2]
    chroms = [tx['chrom']] * len(starts)
    return(list(zip(chroms, ends, starts)))

parser = argparse.ArgumentParser()
parser.add_argument(dest='samfile')
parser.add_argument(dest='outbam_file')
parser.add_argument('--splice_juncs', dest='tx_info', default='',
                    help='''Reference file containing transcripts and their respective
                    splice junctions. Implies that contigs are being filtered against the
                    genome, otherwise transcriptome is assumed.''')
parser.add_argument('--groupings', dest='groupings', default='',
                    help='''Fasta file containing transcriptome that GMAP was aligned to;
                    requires gene symbols to be in the header lines.
                    Used to match transcriptome mappings of novel contigs to genes.
                    When this option is used, an all.groupings file is created instead
                    of an interesting_contigs.txt file.''')
#TODO: get rid of the gen24.info annotation file, as it is redundant. You can get this info from the gtf.

args        = parser.parse_args()
samfile     = args.samfile
outbam_file = args.outbam_file
tx_info     = args.tx_info
txome_fasta    = args.groupings

outbam_file_unsort = '%s_unsorted.bam' % os.path.splitext(outbam_file)[0]
groupings = []

if tx_info != '':
    genref = pd.read_csv(tx_info, sep='\t')
    juncs = genref.apply(lambda tx: get_juncs(tx), axis=1)
    juncs = [(str(c), int(s), int(e)) for jv in juncs.values for c, s, e in jv] # flatten juncs list

lookup = []
if txome_fasta != '':
    for record in SeqIO.parse(txome_fasta, 'fasta'):
        tx_id = record.id
        gname = re.search('gene_symbol:([A-Za-z0-9\.\_\-]+)', record.description).group(1)
        lookup.append((tx_id, gname))
    lookup = pd.DataFrame(lookup, columns=['tx_id', 'gene_name'])

sam = pysam.AlignmentFile(samfile, 'rc')
outbam = pysam.AlignmentFile(outbam_file_unsort, 'wb', template=sam)

# write novel contigs to bam file
int_contigs = {}
for read in sam.fetch():
    if read.reference_id < 0:
        # unmapped contig
        if read.query_name in int_contigs.keys():
            int_contigs[read.query_name] = np.append(int_contigs[read.query_name], None)
        else:
            int_contigs[read.query_name] = None
        outbam.write(read)
        continue

    has_gaps = any([op in gaps.values() and val >= gap_min for op, val in read.cigar])
    has_clips = any([op in clips.values() and val >= clip_min for op, val in read.cigar])

    unknown_juncs = False
    if tx_info != '':
        # check junctions
        starts, ends = zip(*read.blocks)
        chroms = [read.reference_name] * (len(starts)-1)
        tx_juncs = list(zip(chroms, ends[:-1], starts[1:]))
        unknown_juncs = any([txj not in juncs for txj in tx_juncs])

    if has_gaps or has_clips or unknown_juncs:
        if read.query_name in int_contigs.keys():
            int_contigs[read.query_name] = np.append(int_contigs[read.query_name], read.reference_name)
        else:
            int_contigs[read.query_name] = np.array([read.reference_name])
        outbam.write(read)
sam.close()
outbam.close()

if len(lookup) > 0:
    groupings = []
    all_gns = []
    for contig in int_contigs:
        enst = int_contigs[contig]
        if not enst: continue
        genes = '|'.join(lookup[lookup.tx_id.isin(enst)].gene_name.values)
        all_gns.append(genes)
        groupings.append([contig, genes])
    #TODO: finish code
else:
    # write interesting contigs list to file
    int_contigs = int_contigs.keys()
    outdir = os.path.dirname(outbam_file)
    outdir = '.' if outdir == '' else outdir
    with open('%s/interesting_contigs.txt' % outdir, 'w') as fout:
        for contig in int_contigs:
            fout.write('%s\n' % contig)

pysam.sort('-o', outbam_file, outbam_file_unsort)
pysam.index(outbam_file)
os.remove(outbam_file_unsort)
