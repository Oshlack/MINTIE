#########################################################
# Author: Marek Cmero
# Take an samfile of aligned contigs, and filter out any
# contigs that match annotated transcript characteristics
# i.e. - genomic gaps <7bp
#      - no soft or hard clips >30bp
#      - all junctions are known (in annotation)
# Also groups transcripts to genes (optionally)
#########################################################

import argparse
import pandas as pd
import numpy as np
import pysam
import os
import re
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument(dest='samfile')
parser.add_argument(dest='outbam_file')
parser.add_argument('--splice_juncs', dest='tx_info', default='',
                    help='''Reference file containing transcripts and their respective
                    splice junctions. Implies that contigs are being filtered against the
                    genome, otherwise transcriptome is assumed.''')
parser.add_argument('--groupings', dest='groupings', default='',
                    help='''Transcriptome reference that GMAP was aligned to.
                    Requires gene symbols to be in the header lines.
                    Used to match transcriptome mappings of novel contigs to genes.
                    When this option is used, an all.groupings file is created instead
                    of an interesting_contigs.txt file.''')
parser.add_argument('--annotate', dest='annotate', action='store_true',
                    help='''Annotate all contigs in supplied SAM file''')

args        = parser.parse_args()
samfile     = args.samfile
outbam_file = args.outbam_file
tx_info     = args.tx_info
txome_fasta = args.groupings
annotate    = args.annotate

# cutoff parameters
gap_min = 7
clip_min = 30
match_min = 30
match_perc_min = 0.3

# CIGAR specification codes
gaps = {1: 'insertion', 2: 'deletion', 6: 'silent_deletion'}
clips = {4: 'soft', 5: 'hard'}
outbam_file_unsort = '%s_unsorted.bam' % os.path.splitext(outbam_file)[0]
groupings = []

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

def annotate_contig(read, tx_juncs):
    '''
    return the putative cryptic variant
    and the affected genomic positions
    '''
    gap_idxs = [idx for idx, gap in enumerate(read.cigar) if gap[0] in gaps and gap[1] >= gap_min]
    clip_idxs = [idx for idx, clip in enumerate(read.cigar) if clip[0] in clips and clip[1] >= gap_min]
    contig_juncs = [] if tx_juncs is None else [txj for txj in tx_juncs if txj not in juncs]

    annot = []
    match_idxs = [idx for idx,cig in enumerate(read.cigar) if cig[0] == 0]
    blocks = [b for b in zip(match_idxs, read.blocks)]
    chrom1, chrom2 = read.reference_name, read.reference_name

    if len(gap_idxs) > 0:
        for gap_idx in gap_idxs:
            cigar = read.cigar[gap_idx]
            gtype, size = gaps[cigar[0]], int(cigar[1])

            block_idx = 0 if gap_idx == 0 else np.max(np.where(np.array([b[0] for b in blocks])<gap_idx)[0])
            block = blocks[block_idx][1]
            pos1 = int(block[1])
            pos2 = pos1 if gap_idx == 0 else int(read.blocks[block_idx][1])

            annot.append([read.query_name, gtype, chrom1, pos1, chrom2, pos2, size])
            print('%d %s at pos %s:%d-%d (cigar string = %s)' % (size, gtype, chrom1, pos1, pos2, read.cigarstring))

    if len(clip_idxs) > 0:
        for clip_idx in clip_idxs:
            cigar = read.cigar[clip_idx]
            gtype = 'fusion' if cigar[0]==5 else 'soft-clip'

            block_idx = 0 if clip_idx == 0 else np.max(np.where(np.array([b[0] for b in blocks])<clip_idx)[0])
            block = read.blocks[block_idx]
            pos1, size = block[1], read.cigar[clip_idx][1]

            annot.append([read.query_name, gtype, chrom1, pos1, 'NA', 0, size])
            print('%d bp %s at pos %s:%d (cigar string = %s)' % (size, gtype, chrom1, pos1, read.cigarstring))

    if len(contig_juncs) > 0:
        for junc in contig_juncs:
            pos1, pos2 = int(junc[1]), int(junc[2])
            annot.append([read.query_name, 'novel junction', chrom1, pos1, chrom2, pos2, pos2 - pos1])
            print('novel junction at pos %s:%s-%s (cigar string = %s)' % (chrom1, junc[1], junc[2], read.cigarstring))

    return(annot)

if tx_info != '':
    # aligning against genome
    genref = pd.read_csv(tx_info, sep='\t')
    juncs = genref.apply(lambda tx: get_juncs(tx), axis=1)
    juncs = [(str(c), int(s), int(e)) for jv in juncs.values for c, s, e in jv] # flatten juncs list
else:
    # this means we are analysing against transcriptome (not genome)
    gaps = {'insertion': 1, 'deletion': 2, 'skipped': 3, 'silent_deletion': 6} # also consider skipped regions as gaps

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
int_contigs = []
for read in sam.fetch():
    if read.reference_id < 0 or read.mapping_quality == 0:
        # skip unmapped or 0 MAPQ contigs
        continue

    # only consider the contig if at least match_min bases align
    # to reference and at least match_perc_min of the read aligns
    rlen = read.reference_length
    qlen = float(read.query_length)
    if (rlen < match_min) or (rlen / qlen) < match_perc_min:
        continue

    has_gaps = any([op in gaps.keys() and val >= gap_min for op, val in read.cigar])
    has_clips = any([op in clips.keys() and val >= clip_min for op, val in read.cigar])

    tx_juncs = None
    unknown_juncs = False
    if tx_info != '':
        # check junctions
        starts, ends = zip(*read.blocks)
        chroms = [read.reference_name] * (len(starts)-1)
        tx_juncs = list(zip(chroms, ends[:-1], starts[1:]))
        unknown_juncs = any([txj not in juncs for txj in tx_juncs])

    if has_gaps or has_clips or unknown_juncs:
        if annotate:
            annotation = annotate_contig(read, tx_juncs)
            int_contigs.extend(annotation)
        else:
            int_contigs.append([read.query_name, 'novel_contig', read.reference_name])
        outbam.write(read)

sam.close()
outbam.close()

outdir = os.path.dirname(outbam_file)
outdir = '.' if outdir == '' else outdir

if len(lookup) > 0:
    config_dict = {}
    for annot in int_contigs:
        contig = annot[0]
        if contig in contig_dict:
           contig_dict[contig] = np.append(contig_dict[contig], annot[2])
        else:
           contig_dict[contig] = np.array(annot[2])

    groupings = pd.DataFrame()
    all_gns = pd.DataFrame()

    for contig in contig_dict:
        enst = contig_dict[contig]
        if enst is None: continue
        genes = np.unique(lookup[lookup.tx_id.isin(enst)].gene_name.values)
        fus_genes = '|'.join(genes)
        for gn in genes:
            all_gns = all_gns.append([[gn, fus_genes]])
        groupings = groupings.append([[contig, fus_genes]])

    # sort contigs alphanumerically
    tmp = [int(x.split('_')[1]) for x in groupings[0].values]
    groupings = groupings.reset_index()
    groupings = groupings.drop(['index'], axis=1)
    groupings = groupings.loc[sorted(range(len(tmp)), key=lambda k: tmp[k])]

    groupings = pd.concat([groupings, all_gns.sort_values(by=0)])
    groupings = groupings.drop_duplicates()
    groupings.to_csv('%s/all.groupings' % outdir, sep='\t', header=False, index=False)
else:
    # write interesting contigs list to file
    write_header = False
    if annotate:
        int_contigs = pd.DataFrame(int_contigs, columns=['contig', 'variant', 'chrom1', 'pos1', 'chrom2', 'pos2', 'size'])
        write_header = True
    else:
        int_contigs = np.unique(np.array([c[0] for c in contigs]))
        int_contigs = pd.DataFrame(list(int_contigs))
    int_contigs.to_csv('%s/interesting_contigs.txt' % outdir, sep='\t', header=write_header, index=False)

pysam.sort('-o', outbam_file, outbam_file_unsort)
pysam.index(outbam_file)
os.remove(outbam_file_unsort)
