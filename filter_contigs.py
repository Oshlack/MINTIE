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
parser.add_argument(dest='samfile',
                    help='''SAM or BAM format file containing contig alignments''')
parser.add_argument(dest='outbam_file',
                    help='''BAM file to write contigs which pass filtering''')
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
parser.add_argument('--annotate', dest='annotate', default='',
                    help='''Differential splicing results file. This option will annotate
                    all contigs in the supplied SAM file, and write an annotation file
                    containing the novel variant(s) per contig.''')

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

sam = pysam.AlignmentFile(samfile, 'rc')
outbam = pysam.AlignmentFile(outbam_file_unsort, 'wb', template=sam)

ds_output = None
if annotate != '':
    ds_output = pd.read_csv(annotate, sep='\t')

def get_juncs(tx):
    '''
    return list of junctions in form
    [(chr, start, end)] from transcript info file
    note that exon *ends* become junction *starts*
    '''
    starts = tx['exonStarts'].split(',')[1:]
    ends = tx['exonEnds'].split(',')[:-1]
    chroms = [tx['chrom']] * len(starts)
    return(list(zip(chroms, ends, starts)))

def annotate_contig(read, tx_juncs):
    '''
    return the putative cryptic variant
    and the affected genomic positions
    '''
    gap_idxs = [idx for idx, gap in enumerate(read.cigar) if gap[0] in gaps and gap[1] >= gap_min]
    clip_idxs = [idx for idx, clip in enumerate(read.cigar) if clip[0] in clips and clip[1] >= gap_min]
    #contig_juncs = [] if tx_juncs is None else [txj for txj in tx_juncs if txj not in juncs]

    annot = []
    match_idxs = [idx for idx,cig in enumerate(read.cigar) if cig[0] == 0]
    blocks = [b for b in zip(match_idxs, read.get_blocks())]
    chrom1, chrom2 = read.reference_name, read.reference_name

    if len(gap_idxs) > 0:
        for gap_idx in gap_idxs:
            cigar = read.cigar[gap_idx]
            gtype, size = gaps[cigar[0]], int(cigar[1])

            block_idx = 0 if gap_idx == 0 else np.max(np.where(np.array([b[0] for b in blocks])<gap_idx)[0])
            block = blocks[block_idx][1]
            pos1 = int(block[1])
            pos2 = pos1
            if gap_idx != 0:
                pos2 = int(read.blocks[block_idx+1][0]) if block_idx+1 < len(blocks) else pos1 + size

            # position of variant on contig
            cpos1 = sum([v for c,v in read.cigar[:gap_idx]])
            csize = size if read.cigar[gap_idx][0] == 1 else 0 # only insertions affect contig pos
            cpos2 = cpos1 + csize

            var_seq = ''
            if csize > 0:
                seq_pos1 = sum([v for c,v in read.cigar[:gap_idx] if c!=5])
                seq_pos2 = seq_pos1 + csize
                var_seq = read.query_sequence[seq_pos1:seq_pos2]

            annot.append([read.query_name, gtype, chrom1, pos1, chrom2, pos2, size, cpos1, cpos2, csize, var_seq])
            print('%d %s at pos %s:%d-%d (cigar string = %s)' % (size, gtype, chrom1, pos1, pos2, read.cigarstring))

    if len(clip_idxs) > 0:
        for clip_idx in clip_idxs:
            cigar = read.cigar[clip_idx]
            gtype = 'fusion' if cigar[0]==5 else 'soft-clip'

            block_idx = 0 if clip_idx == 0 else np.max(np.where(np.array([b[0] for b in blocks])<clip_idx)[0])
            block = read.get_blocks()[block_idx]
            pos1, size = block[1], read.cigar[clip_idx][1]

            # position of variant on contig
            cpos1 = sum([v for c,v in read.cigar[:clip_idx]])
            csize = 0 if gtype == 'fusion' else cigar[1]
            cpos2 = cpos1 if gtype == 'fusion' else cpos1 + csize

            var_seq = ''
            if csize > 0:
                seq_pos1 = sum([v for c,v in read.cigar[:clip_idx] if c!=5])
                seq_pos2 = seq_pos1 + csize
                var_seq = read.query_sequence[seq_pos1:seq_pos2]

            annot.append([read.query_name, gtype, chrom1, pos1, chrom1, pos1, size, cpos1, cpos2, csize, var_seq])
            print('%d bp %s at pos %s:%d (cigar string = %s)' % (size, gtype, chrom1, pos1, read.cigarstring))

    if len(tx_juncs) > 0:
        for junc in tx_juncs:
            # don't consider if the junction is a deletion, insertion or clipped sequence
            pos1, pos2 = int(junc[1]), int(junc[2])
            junc_idx = [idx for idx, block in blocks if block[1] == pos1][0]
            junc_type = read.cigar[junc_idx+1][0]
            if junc_type in gaps or junc_type in clips:
                continue

            # position of variant on contig
            cpos1 = sum([v for c,v in read.cigar[:(junc_idx+1)]])
            cpos2 = cpos1

            #TODO: extract intra/intergenic sequence from novel junction
            annot.append([read.query_name, 'novel junction', chrom1, pos1, chrom2, pos2, pos2 - pos1, cpos1, cpos2, 0, ''])
            print('novel junction at pos %s:%s-%s (cigar string = %s)' % (chrom1, junc[1], junc[2], read.cigarstring))

    return(annot)

def nice_sort(ids):
    '''
    sorts strings with characters and integers intuitively
    '''
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(ids, key = alphanum_key)

def pair_fusions(novel_contigs):
    fusions = novel_contigs[novel_contigs.variant == 'fusion']
    fusion_contigs = np.unique(fusions.contig).copy()
    if len(fusion_contigs) == 0:
        return novel_contigs

    fusions_new = pd.DataFrame()
    for contig in fusion_contigs:
        tmp = fusions[fusions.contig==contig].reset_index(drop=True)
        if len(tmp) != 2:
            continue
        chroms, pos = nice_sort(tmp.chrom1.values), np.sort(tmp.genome_pos1.values)
        if tmp.chrom1[0] != tmp.chrom1[1]:
            genome_pos1 = tmp[chroms[0] == tmp.chrom1].genome_pos1.values[0]
            genome_pos2 = tmp[chroms[1] == tmp.chrom1].genome_pos1.values[0]
            pos = [genome_pos1, genome_pos2]
        tmp.loc[0, 'chrom1'] = chroms[0]
        tmp.loc[0, 'chrom2'] = chroms[1]
        tmp.loc[0, 'genome_pos1'] = pos[0]
        tmp.loc[0, 'genome_pos2'] = pos[1]
        fusions_new = pd.concat([fusions_new, tmp.loc[0]], axis=1)

    fusions_new = fusions_new.transpose()

    # set size to correspond to rearrangement proximity
    fusions_new['size'] = None
    intra_chrom = fusions_new.chrom1 == fusions_new.chrom2
    fusions_new.loc[intra_chrom, 'size'] = fusions_new[intra_chrom].genome_pos2.values - fusions_new[intra_chrom].genome_pos1.values

    novel_contigs = pd.concat([novel_contigs[novel_contigs.variant != 'fusion'], fusions_new])
    return(novel_contigs.reset_index(drop=True))

if tx_info != '':
    print('Generating lookup for known splice junctions...')
    # aligning against genome
    genref = pd.read_csv(tx_info, sep='\t')
    juncs = genref.apply(lambda tx: get_juncs(tx), axis=1)

    #juncs = [(str(c), int(s), int(e)) for jv in juncs.values for c, s, e in jv] # flatten juncs list
    juncs = ['%s:%s-%s' % (c, s, e) for jv in juncs.values for c, s, e in jv] # flatten juncs list
    junc_dic = {}
    for junc in juncs:
        junc_dic[junc] = True
    juncs = junc_dic
    print('Finished generating lookup')
else:
    # this means we are analysing against transcriptome (not genome)
    gaps = {1: 'insertion', 2: 'deletion', 3: 'skipped',  6: 'silent_deletion'} # consider skipped regins as gaps

lookup = []
if txome_fasta != '':
    for record in SeqIO.parse(txome_fasta, 'fasta'):
        tx_id = record.id
        gname = re.search('gene_symbol:([A-Za-z0-9\.\_\-]+)', record.description).group(1)
        lookup.append((tx_id, gname))
    lookup = pd.DataFrame(lookup, columns=['tx_id', 'gene_name'])

# write novel contigs to bam file
print('Checking contigs for non-reference content...')
novel_contigs = []
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

    has_gaps = any([op in gaps and val >= gap_min for op, val in read.cigar])
    has_clips = any([op in clips and val >= clip_min for op, val in read.cigar])

    tx_juncs = []
    unknown_juncs = []
    if tx_info != '':
        # check junctions
        starts, ends = zip(*read.get_blocks())
        chroms = [read.reference_name] * (len(starts)-1)
        tx_juncs = list(zip(chroms, ends[:-1], starts[1:]))
        tx_juncs = [junc for junc in tx_juncs if (junc[2] - junc[1]) > gap_min]
        unknown_juncs = ['%s:%s-%s' % (c, s, e) not in juncs for c, s, e in tx_juncs]

    if has_gaps or has_clips or any(unknown_juncs):
        if annotate:
            novel_juncs = [list(x) for x in np.array(tx_juncs)[unknown_juncs]]
            annotation = annotate_contig(read, novel_juncs)
            novel_contigs.extend(annotation)
        else:
            novel_contigs.append([read.query_name, 'novel_contig', read.reference_name])
        outbam.write(read)

sam.close()
outbam.close()

outdir = os.path.dirname(outbam_file)
outdir = '.' if outdir == '' else outdir

print('Writing results...')
if len(lookup) > 0:
    contig_dict = {}
    for annot in novel_contigs:
        contig = annot[0]
        if contig in contig_dict:
           contig_dict[contig] = contig_dict[contig] + [annot[2]]
        else:
           contig_dict[contig] = [annot[2]]

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
    annot = ''
    if annotate != '':
        cols = ['contig', 'variant', 'chrom1', 'genome_pos1', 'chrom2', \
                'genome_pos2', 'size', 'contig_pos1', 'contig_pos2', \
                'contig_varsize', 'variant_seq']
        novel_contigs = pd.DataFrame(novel_contigs, columns=cols)
        novel_contigs = pair_fusions(novel_contigs)
        novel_contigs = novel_contigs.merge(ds_output, left_on='contig', right_on='transcript', how='inner')
        output_cols = ['gene', 'contig', 'variant', 'chrom1', 'genome_pos1',
                       'chrom2', 'genome_pos2', 'size', 'contig_pos1',
                       'contig_pos2', 'contig_varsize', 'variant_seq',
                       'ec_names', 'contigs', 'padj', 'gene.FDR',
                       'UniqueCount', 'AmbigCount', 'ambig_ratio']
        novel_contigs = novel_contigs[output_cols].drop_duplicates()
        novel_contigs = novel_contigs.sort_values(by=['padj'])
        write_header = True
        annot = '_annotated'
    else:
        novel_contigs = np.unique(np.array([c[0] for c in novel_contigs]))
        novel_contigs = pd.DataFrame(list(novel_contigs))
    novel_contigs.to_csv('%s/novel_contigs%s.txt' % (outdir, annot), sep='\t', header=write_header, index=False)

pysam.sort('-o', outbam_file, outbam_file_unsort)
pysam.index(outbam_file)
os.remove(outbam_file_unsort)
