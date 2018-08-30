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
parser.add_argument('--tx_align', dest='tx_align', default='',
                    help='''Alignment of novel contigs to transcriptome. Used for more detailed
                    annotation.''')

args        = parser.parse_args()
samfile     = args.samfile
outbam_file = args.outbam_file
tx_info     = args.tx_info
txome_fasta = args.groupings
annotate    = args.annotate
tx_align    = args.tx_align

# cutoff parameters
gap_min = 7
clip_min = 30
match_min = 30
match_perc_min = 0.3

# CIGAR specification codes
gaps = {1: 'insertion', 2: 'deletion', 6: 'silent_deletion'}
clips = {4: 'soft', 5: 'hard'}
ref_only_gaps = {2: 'deletion', 3: 'skipped', 5: 'hard-clip'}
outbam_file_unsort = '%s_unsorted.bam' % os.path.splitext(outbam_file)[0]
groupings = []

sam = pysam.AlignmentFile(samfile, 'rc')
outbam = pysam.AlignmentFile(outbam_file_unsort, 'wb', template=sam)
tx_bam = None
if tx_align != '':
    tx_bam = pysam.AlignmentFile(tx_align, 'rc')
    tx_idx = pysam.IndexedReads(tx_bam)
    tx_idx.build()

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
    contig_size = sum([v for c,v in read.cigar if c in [0, 1, 4, 5]]) # count if match, insertion or clip

    annot, var_seq = [], ''
    match_idxs = [idx for idx,cig in enumerate(read.cigar) if cig[0] == 0]
    blocks = [b for b in zip(match_idxs, read.get_blocks())]
    chrom1, chrom2 = read.reference_name, read.reference_name
    strand = '-' if read.is_reverse else '+'
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
                seq_pos1 = sum([v for c,v in read.cigar[:gap_idx] if c not in ref_only_gaps])
                seq_pos2 = seq_pos1 + csize
                var_seq = read.query_sequence[seq_pos1:seq_pos2]

            annot.append([read.query_name, gtype, chrom1, pos1, chrom2, pos2, contig_size, size, cpos1, cpos2, strand, csize, var_seq])
            print('%d %s at pos %s:%d-%d (cigar string = %s)' % (size, gtype, chrom1, pos1, pos2, read.cigarstring))

    if len(clip_idxs) > 0:
        for clip_idx in clip_idxs:
            cigar = read.cigar[clip_idx]
            gtype = 'fusion' if cigar[0]==5 else 'soft-clip'

            block_idx = 0 if clip_idx == 0 else np.max(np.where(np.array([b[0] for b in blocks])<clip_idx)[0])
            block = read.get_blocks()[block_idx]
            size = read.cigar[clip_idx][1]
            pos1 = block[1] if clip_idx > 0 else block[0] # pick coord based on which side clip is on

            # position of variant on contig
            cpos1 = sum([v for c,v in read.cigar[:clip_idx]])
            cpos1 = contig_size - size if ((read.is_reverse and clip_idx == 0) or (not read.is_reverse and clip_idx > 0)) else size
            csize = 0 if gtype == 'fusion' else cigar[1]
            cpos2 = cpos1 if gtype == 'fusion' else cpos1 + csize

            var_seq = ''
            if csize > 0:
                seq_pos1 = sum([v for c,v in read.cigar[:clip_idx] if c not in ref_only_gaps])
                seq_pos2 = seq_pos1 + csize
                var_seq = read.query_sequence[seq_pos1:seq_pos2]

            if gtype == 'fusion' and tx_bam:
                tx_reads = [tx_read for tx_read in tx_idx.find(read.query_name)]
                if len(tx_reads) == 1:
                    tx_read = tx_reads[0]
                    is_softclip = [o == 4 for o,v in tx_read.cigar]

                    if any(is_softclip):
                        sc_start = is_softclip[0]
                        if is_softclip[0] and is_softclip[-1]:
                            sc_start = tx_read.cigar[0][1] > tx_read.cigar[-1][1]
                            print('WARNING: contig %s is soft-clipped at both ends. Picking the longest soft-clip for annotation.' % read.query_name)

                        var_seq = str(tx_read.query_sequence)
                        #sc_size = tx_read.cigar[0][1] if sc_start else tx_read.cigar[-1][1]
                        sc_size = read.cigar[clip_idx][1]
                        # ^ actually the hard clip len from genome alignment, pick this as the sc_size
                        # (even though it may differ from the txome alignment), otherwise we don't know
                        # the genomic pos where novel seq is inserted (need this to build supertranscript)
                        var_seq = var_seq[:sc_size] if sc_start else var_seq[-sc_size:]
                        cpos1 = 0 if sc_start else len(tx_read.query_sequence) - sc_size
                        cpos2 = cpos1 + sc_size

            annot.append([read.query_name, gtype, chrom1, pos1, chrom1, pos1, contig_size, size, cpos1, cpos2, strand, csize, var_seq])
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

            if not tx_bam:
                annot.append([read.query_name, 'novel junction', chrom1, pos1, chrom2, pos2,
                              contig_size, pos2 - pos1, cpos1, cpos2, strand, 0, ''])
                continue

            tx_reads = [tx_read for tx_read in tx_idx.find(read.query_name)]
            if len(tx_reads) > 2:
                print('WARNING: cannot annotate contig %s because it maps to >2 transcripts.' % read.query_name)
                annot.append([read.query_name, 'novel junction; multimap', chrom1, pos1, chrom2, pos2,
                              contig_size, pos2 - pos1, cpos1, cpos2, strand, 0, ''])
                continue

            for tx_read in tx_reads:
                is_softclip = [o == 4 for o,v in tx_read.cigar]
                if not any(is_softclip):
                    continue
                sc_start = is_softclip[0]
                if is_softclip[0] and is_softclip[-1]:
                    sc_start = tx_read.cigar[0][1] > tx_read.cigar[-1][1]
                    #TODO: include both ends?
                    print('WARNING: contig %s is soft-clipped at both ends. Picking the longest soft-clip for annotation.' % read.query_name)
                var_seq = str(tx_read.query_sequence)
                sc_size = tx_read.cigar[0][1] if sc_start else tx_read.cigar[-1][1]
                var_seq = var_seq[:sc_size] if sc_start else var_seq[-sc_size:]
                cpos1 = 0 if sc_start else len(tx_read.query_sequence) - sc_size
                cpos2 = cpos1 + sc_size

            csize = len(var_seq)
            variant_type = 'novel splice junction' if csize == 0 else 'junction with novel seq'
            if not csize > 0:
                junc_loc_left = '%s:%d' % (chrom1, pos1) in locs
                junc_loc_right = '%s:%d' % (chrom2, pos2) in locs
                if not junc_loc_left and not junc_loc_right:
                    variant_type = 'non-exonic junction'
                elif not (junc_loc_left and junc_loc_right):
                    variant_type = 'exonic-to-non-exonic junction'

            annot.append([read.query_name, variant_type, chrom1, pos1, chrom2, pos2, contig_size, pos2 - pos1, cpos1, cpos2, strand, csize, var_seq])
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
    junc_info = genref.apply(lambda tx: get_juncs(tx), axis=1)

    #juncs = [(str(c), int(s), int(e)) for jv in juncs.values for c, s, e in jv] # flatten juncs list
    juncs = ['%s:%s-%s' % (c, s, e) for jv in junc_info.values for c, s, e in jv] # flatten juncs list
    junc_dic = {}
    for junc in juncs:
        junc_dic[junc] = True
    juncs = junc_dic

    # create a list of junction start/ends for more detailed annotation
    locs = [['%s:%s' % (c, s), '%s:%s' % (c, e)] for jv in junc_info.values for c, s, e in jv]
    locs = [l for loc in locs for l in loc] #unlist
    loc_dic = {}
    for loc in locs:
        loc_dic[loc] = True
    locs = loc_dic
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
    if read.reference_id < 0:
        # skip unmapped contigs
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
if tx_align != '' :
    tx_bam.close()

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
                'genome_pos2', 'contig_size', 'genome_varsize', 'contig_pos1', 'contig_pos2', \
                'contig_align_strand', 'contig_varsize', 'variant_seq']
        novel_contigs = pd.DataFrame(novel_contigs, columns=cols)
        novel_contigs = pair_fusions(novel_contigs)
        novel_contigs = novel_contigs.merge(ds_output, left_on='contig', right_on='transcript', how='inner')
        sample = os.path.dirname(annotate).split('/')[-1].split('_')[0]
        novel_contigs['sample'] = sample
        output_cols = ['gene', 'contig', 'variant', 'chrom1', 'genome_pos1',
                       'chrom2', 'genome_pos2', 'contig_size', 'genome_varsize',
                       'contig_pos1', 'contig_pos2', 'contig_align_strand',
                       'contig_varsize', 'variant_seq', 'ec_names',
                       'contigs', 'FDR', 'UniqueCount',
                       'AmbigCount', 'ambig_ratio', 'sample']
        novel_contigs = novel_contigs[output_cols].drop_duplicates()
        novel_contigs = novel_contigs.sort_values(by=['FDR'])
        write_header = True
        annot = '_annotated'
    else:
        novel_contigs = np.unique(np.array([c[0] for c in novel_contigs]))
        novel_contigs = pd.DataFrame(list(novel_contigs))
    novel_contigs.to_csv('%s/novel_contigs%s.txt' % (outdir, annot), sep='\t', header=write_header, index=False)

pysam.sort('-o', outbam_file, outbam_file_unsort)
pysam.index(outbam_file)
os.remove(outbam_file_unsort)
