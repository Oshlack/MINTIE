'''
Module      : annotate_contigs
Description : Perform contig annotation.
Copyright   : (c) Marek Cmero, Sep 2018
License     : TBD
Maintainer  : MAREK.CMERO@MCRI.EDU.AU
Portability : POSIX
Take an samfile of aligned contigs, filter out
any contigs matching the reference annotation,
and annotate variant(s) contained in each contig
with novel component(s).
'''

import pandas as pd
import numpy as np
import pysam
import os
import pickle
import re
import logging
import sys
import string
from argparse import ArgumentParser
from intervaltree import Interval, IntervalTree
from cv_vcf import CrypticVariant, VCF
from utils import cached, init_logging, exit_with_error
import ipdb

pd.set_option("mode.chained_assignment", None)

EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
PROGRAM_NAME = "annotate_contigs"

# CUT-OFF parameters
GAP_MIN = 7
CLIP_MIN = 30
MATCH_MIN = 30
MATCH_PERC_MIN = 0.3

# VCF parameters
INFO = ["CID", "ECN", "CLEN", "CPOS", "CSTRAND", "CCIGAR", "VSIZE", "CVSIZE", "CVTYPE", "GENES", "PARID", "PVAL", "CVQ"]
FORMAT = ["GT", "ECC", "AI"]

# CIGAR specification codes
CIGAR = {'match': 0,
         'insertion': 1,
         'deletion': 2,
         'skipped': 3,
         'soft-clip': 4,
         'hard-clip': 5,
         'silent_deletion': 6}
GAPS = [CIGAR[c] for c in ['insertion', 'deletion', 'silent_deletion']]
CLIPS = [CIGAR[c] for c in ['soft-clip', 'hard-clip']]
# any cigar criteria that is >0 bp on an aligned contig
AFFECT_CONTIG = [CIGAR[c] for c in ['match', 'insertion', 'soft-clip', 'hard-clip']]
# any cigar criteria that is >0 bp on the reference genome
AFFECT_REF = [CIGAR[c] for c in ['match', 'deletion']]

# read processing record
record = {}


def parse_args():
    '''
    Parse command line arguments.
    Returns Options object with command line argument values as attributes.
    Will exit the program on a command line error.
    '''
    description = 'Annotate contigs'
    parser = ArgumentParser(description=description)
    parser.add_argument('--log',
                        metavar='LOG_FILE',
                        type=str,
                        help='record program progress in LOG_FILE')
    parser.add_argument(dest='sample',
                        metavar='SAMPLE',
                        type=str,
                        help='''Sample name''')
    parser.add_argument(dest='bam_file',
                        metavar='BAM_FILE',
                        type=str,
                        help='''SAM or BAM format file containing contig alignments''')
    parser.add_argument(dest='output_bam',
                        metavar='OUTPUT_BAM',
                        type=str,
                        help='''BAM file to write contigs which pass filtering''')
    parser.add_argument(dest='contig_info_file',
                        metavar='CONTIG_INFO_FILE',
                        type=str,
                        help='''Contig info output file''')
    parser.add_argument(dest='junc_file',
                        metavar='JUNC_FILE',
                        type=str,
                        help='''Reference file containing transcripts and their
                        respective splice junctions.''')
    parser.add_argument(dest='tx_ref_file',
                        type=str,
                        metavar='TX_REF_FILE',
                        help='''Transcriptiome GTF reference file.''')
    return parser.parse_args()

def get_gene(attribute):
    re_gene = re.search('gene_name "([\w\-\.\/]+)"', attribute)
    gene = re_gene.group(1) if re_gene else ''
    return gene

@cached('gene_lookup_cache.pickle')
def get_gene_lookup(tx_ref_file):
    '''
    Generate start/end coordinate reference
    for genes and output as an interval tree
    dictionary. Also output dataframe containing
    chromosome, start and ends for all exons.
    '''
    ref_trees, ex_ref_out = None, None
    if tx_ref_file != '':
        logging.info('Generating lookup for genes...')

        #TODO: standardise with make_supertranscript for gtf handling
        tx_ref = pd.read_csv(tx_ref_file, comment='#', sep='\t', header=None)
        tx_ref['gene'] = tx_ref[8].apply(lambda x: get_gene(x))
        aggregator = {3: lambda x: min(x),
                      4: lambda x: max(x)}
        gn_ref = tx_ref.groupby([0, 'gene'], as_index=False, sort=False).agg(aggregator)
        gn_ref = gn_ref[[0, 3, 4, 'gene']]
        gn_ref.columns = ['chrom', 'start', 'end', 'gene']
        gn_ref = gn_ref.drop_duplicates()

        # start/end coordinates for gene matching
        ref_trees = {}
        chroms = np.unique(gn_ref.chrom.values)
        for chrom in chroms:
            chr_ref = gn_ref[gn_ref.chrom == chrom].drop_duplicates()
            ref_tree = IntervalTree()
            for s,e,g in zip(chr_ref['start'].values, chr_ref['end'].values, chr_ref['gene'].values):
                ref_tree.addi(s-1, e, g)
            ref_trees[chrom] = ref_tree

        # merged exon boundaries for block annotation
        ex_ref = tx_ref[tx_ref[2] == 'exon']
        ex_ref_out = pd.DataFrame()
        ex_trees = {}
        for chrom in chroms:
            chr_ref = ex_ref[ex_ref[0] == chrom].drop_duplicates()
            ex_tree = IntervalTree()
            for s,e in zip(chr_ref[3].values, chr_ref[4].values):
                ex_tree.addi(s-1, e)
            ex_tree.merge_overlaps()
            tmp = pd.DataFrame([(chrom, tree[0], tree[1]) for tree in ex_tree],
                               columns=['chrom', 'start', 'end'])
            ex_ref_out = pd.concat([ex_ref_out, tmp], ignore_index=True)
            ex_trees[chrom] = ex_tree
    return ref_trees, ex_trees, ex_ref_out

def get_juncs(tx):
    '''
    Return list of junctions in form
    [(chr, start, end)] from transcript info file
    note that exon *ends* become junction *starts*
    '''
    starts = tx['exonStarts'].split(',')[1:]
    ends = tx['exonEnds'].split(',')[:-1]
    chroms = [tx['chrom']] * len(starts)
    return(list(zip(chroms, ends, starts)))

@cached('juncs_lookup_cache.pickle')
def get_junc_lookup(junc_file):
    '''
    Take junction reference file, generate and
    return two dictionaries, one containing all
    splice junction start and end sites, and one
    containing all coordinates where a junction
    starts or ends.
    '''
    logging.info('Generating lookup for known splice junctions...')

    # aligning against genome
    genref = pd.read_csv(junc_file, sep='\t')
    junc_info = genref.apply(lambda tx: get_juncs(tx), axis=1)

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
    logging.info('Finished generating splice junction lookup')

    return (juncs, locs)

def get_overlapping_genes(read, ref_trees):
    try:
        ref_tree = ref_trees[read.reference_name]
    except KeyError:
       logging.info('WARNING: reference chromosome %s (from read %s) was not found in supplied reference' %
                    (read.reference_name, read.query_name))
       return ''

    blocks = read.get_blocks()
    genes = []
    for block in blocks:
        bstart, bend = block
        gns = ref_tree.overlap(bstart, bend)
        genes.extend([gn[2] for gn in gns if gn[2] != ''])
    genes = np.unique(np.array(genes))

    return('|'.join(genes))

def get_next_letter(last_letter):
    '''
    Convenience function to get next letter in alphabet
    '''
    try:
        next_letter_pos = np.where(np.array(list(string.ascii_letters)) == last_letter)[0][0]+1
        next_letter = list(string.ascii_letters)[next_letter_pos]
        return next_letter
    except IndexError:
        logging.info('Cannot increment letter, non-ascii letter (%s) provided or next letter out of bounds' % last_letter)
        return ''

def get_next_id(qname):
    if qname not in record:
        vid = qname + 'a'
        record[qname] = [vid]
        return vid
    else:
        vid = qname + get_next_letter(record[qname][-1][-1])
        record[qname].append(vid)
        return vid

def annotate_gaps(cv, read, ci_file):
    '''
    Annotate deletions and insertions
    '''
    gap_idxs = [idx for idx, gap in enumerate(read.cigar) if gap[0] in GAPS and gap[1] >= GAP_MIN]
    for gap_idx in gap_idxs:
        cigar = read.cigar[gap_idx]
        cv.vsize = int(cigar[1])

        block_idx = 0 if gap_idx == 0 else np.max(np.where(np.array([b[0] for b in cv.blocks]) < gap_idx)[0])
        block = cv.blocks[block_idx][1]
        cv.pos = int(block[1])

        # position of variant on contig
        cv.cpos = sum([v for c,v in read.cigar[:gap_idx] if c in AFFECT_CONTIG])
        cv.cvsize = cv.vsize if read.cigar[gap_idx][0] == CIGAR['insertion'] else 0 # only insertions affect contig pos

        if cigar[0] == CIGAR['insertion']:
            cv.cvtype = 'INS'
            seq_pos1 = sum([v for c,v in read.cigar[:gap_idx] if c in AFFECT_CONTIG and c != CIGAR['hard-clip']])
            seq_pos2 = seq_pos1 + cv.cvsize
            seq = read.query_sequence[(seq_pos1-1):seq_pos2]
            cv.ref, cv.alt = seq[:1], seq
        else:
            cv.cvtype = 'DEL'
            seq_pos1 = sum([v for c,v in read.cigar[:gap_idx] if c in AFFECT_REF])
            seq_pos2 = seq_pos1 + cv.vsize
            seq = read.get_reference_sequence()[(seq_pos1-1):seq_pos2]
            cv.ref, cv.alt  = seq, seq[:1]

        cv.vid = get_next_id(read.query_name)
        CrypticVariant.write_contig_info(ci_file, cv)
        print(cv.vcf_output())

def annotate_softclips(cv, read, ci_file):
    sc_idxs = [idx for idx, clip in enumerate(read.cigar) if clip[0] == CIGAR['soft-clip'] and clip[1] >= CLIP_MIN]
    for sc_idx in sc_idxs:
        cigar = read.cigar[sc_idx]
        cv.cvtype, cv.vsize, cv.cvsize = 'UN', int(cigar[1]), int(cigar[1])
        sc_left = sc_idx == 0

        block_idx = 0 if sc_idx == 0 else np.max(np.where(np.array([b[0] for b in cv.blocks]) < sc_idx)[0])
        block = cv.blocks[block_idx][1]
        cv.pos = int(block[0]) + 1 if sc_left else int(block[1])

        rcigar = read.cigar[::-1] if cv.strand == '-' else read.cigar
        cv.cpos = sum([v for c,v in rcigar[:sc_idx] if c in AFFECT_CONTIG])

        varseq = read.query_sequence[cv.cpos:(cv.cpos+cv.cvsize)]
        refseq = read.get_reference_sequence()
        cv.ref = refseq[0] if sc_left else refseq[-1]
        cv.alt = '%s]%s:%d]%s' % (varseq, cv.chrom, cv.pos-1, cv.ref) if sc_left else \
                 '%s[%s:%d[%s' % (cv.ref, cv.chrom, cv.pos+1, varseq)

        cv.vid = get_next_id(read.query_name)
        CrypticVariant.write_contig_info(ci_file, cv)
        print(cv.vcf_output())

def get_chr_ref(read, ex_ref):
    '''
    get exon coordinates for a given read's alignment chromosome
    '''
    chr_ref = ex_ref[ex_ref.chrom == read.reference_name]
    if len(chr_ref) == 0:
       logging.info('WARNING: reference chromosome %s (from read %s) was not found in supplied reference' %
                    (read.reference_name, read.query_name))
    return chr_ref

def is_novel_block(block, chr_ref):
    '''
    checks a read's sequence blocks and returns
    false if the block matches (or is contained)
    within a referenced exon, and true otherwise
    '''
    match = chr_ref[np.logical_and(block[0] >= chr_ref.start, block[1] <= chr_ref.end)]
    block_size = block[1] - block[0]
    return len(match) == 0 and block_size > CLIP_MIN

def get_block_sequence(read, block_idx):
    '''
    return query and reference sequence of
    specified block
    '''
    qseq = read.query_sequence
    rseq = read.get_reference_sequence()

    cpos1 = sum([v for c,v in read.cigar[:block_idx] if c in AFFECT_CONTIG and c != CIGAR['hard-clip']])
    cpos2 = sum([v for c,v in read.cigar[:block_idx+1] if c in AFFECT_CONTIG and c != CIGAR['hard-clip']])

    rpos1 = sum([v for c,v in read.cigar[:block_idx] if c in AFFECT_REF])
    rpos2 = sum([v for c,v in read.cigar[:block_idx+1] if c in AFFECT_REF])

    return qseq[cpos1:cpos2], rseq[rpos1:rpos2]

def annotate_block_right(cv, read, cpos, olapping, block, block_idx):
    '''
    extended exon to the *right* of the reference sequence
    '''
    qseq, rseq = get_block_sequence(read, block_idx)
    seq_left_pos = block[1] - max(olapping.end)
    cv.ref, cv.alt = rseq[(-seq_left_pos):], ']' + qseq[(-seq_left_pos):]
    cv.cpos, cv.pos = cpos, max(olapping.end) + 1
    cv.vsize, cv.cvsize = abs(len(cv.alt)-1 - len(cv.ref)), len(cv.alt)-1
    return cv

def annotate_block_left(cv, read, cpos, olapping, block, block_idx):
    '''
    extended exon is to the *left* of the reference sequence
    '''
    qseq, rseq = get_block_sequence(read, block_idx)
    seq_right_pos = min(olapping.start) - block[0]
    cv.ref, cv.alt = rseq[:seq_right_pos], qseq[:seq_right_pos] + '['
    cv.cpos, cv.pos = cpos, min(olapping.start) - len(cv.ref) + 1
    cv.vsize, cv.cvsize = abs(len(cv.alt)-1 - len(cv.ref)), len(cv.alt)-1
    return cv

def annotate_blocks(cv, read, chr_ref, ci_file):
    '''
    Annotate any sequence that is outside of exonic regions
    '''
    cv.parid = '.' # blocks don't have pairs
    novel_blocks = [(idx, block) for idx, block in cv.blocks if is_novel_block(block, chr_ref)]
    for block_idx, block in novel_blocks:
        cpos1 = sum([v for c,v in read.cigar[:block_idx] if c in AFFECT_CONTIG])
        cpos2 = sum([v for c,v in read.cigar[:block_idx+1] if c in AFFECT_CONTIG])

        olapping = chr_ref[np.logical_and(block[0] < chr_ref.start, block[1] > chr_ref.end)]

        # whether sequence block is on left or right side of contig block, respectively
        left = chr_ref[np.logical_and(block[1] > chr_ref.start, block[1] <= chr_ref.end)]
        right = chr_ref[np.logical_and(block[0] >= chr_ref.start, block[0] < chr_ref.end)]

        if len(left) > 0 and len(right) > 0:
            # retained intron
            cv.cvtype = 'RI'
            qseq, rseq = get_block_sequence(read, block_idx)
            seq_right_pos = block[1] - min(left.start)
            seq_left_pos = max(right.end) - block[0]

            cv.pos = block[0] + seq_left_pos + 1
            cv.ref = rseq[seq_left_pos:(-seq_right_pos)]
            cv.alt = ']' + qseq[seq_left_pos:(-seq_right_pos)] + '['
            cv.cpos = cpos1 + seq_left_pos

            cv.vsize, cv.cvsize = abs(len(cv.alt)-2 - len(cv.ref)), len(cv.alt)-2
            cv.vid = get_next_id(read.query_name)
        elif len(olapping) > 0:
            # annotate left side
            cv.cvtype = 'EE'
            cv = annotate_block_left(cv, read, cpos2, olapping, block, block_idx)
            cv.vid = get_next_id(read.query_name)

            print(cv.vcf_output())
            CrypticVariant.write_contig_info(ci_file, cv)

            # annotate right side
            cv = annotate_block_right(cv, read, cpos1, olapping, block, block_idx)
            cv.vid = get_next_id(read.query_name)
        elif len(left) > 0:
            # annotate left side
            cv.cvtype = 'EE'
            cv = annotate_block_left(cv, read, cpos2, left, block, block_idx)
            cv.vid = get_next_id(read.query_name)
        elif len(right) > 0:
            # annotate right side
            cv.cvtype = 'EE'
            cv = annotate_block_right(cv, read, cpos1, right, block, block_idx)
            cv.vid = get_next_id(read.query_name)
        else:
            # block does not cross any annotation
            qseq, rseq = get_block_sequence(read, block_idx)
            cv.ref, cv.alt = rseq, '[' + qseq + ']'
            cv.pos, cv.cvtype = block[0], 'NE'
            cv.cpos = cpos1
            cv.vid = get_next_id(read.query_name)
            cv.vsize, cv.cvsize = abs(len(cv.alt)-2 - len(cv.ref)), len(cv.alt)-2

        print(cv.vcf_output())
        CrypticVariant.write_contig_info(ci_file, cv)

def annotate_fusion(args, read, juncs, bam_idx, ex_ref, ref_trees, outbam):
    try:
        r1, r2 = bam_idx.find(read.query_name)
    except ValueError:
        logging.info('WARNING: found >2 reads matching hard-clipped read %s; cannot process' % read.query_name)
        return

    ci_file = args.contig_info_file
    cv1 = CrypticVariant().from_read(r1)
    cv2 = CrypticVariant().from_read(r2)
    cv1.cvtype, cv2.cvtype = 'FUS', 'FUS'

    cv1.genes = get_overlapping_genes(r1, ref_trees)
    cv2.genes = get_overlapping_genes(r2, ref_trees)
    if cv1.genes == '' and cv2.genes == '':
        # no intersecting gene, this is not an interesting fusion
        logging.info('No gene(s) intersecting candidate fusion contig %s; skipping' % read.query_name)
        record[read.query_name] = []
        return

    hc_idx1 = [idx for idx, clip in enumerate(r1.cigar) if clip[0] == CIGAR['hard-clip']][0]
    hc_left1 = hc_idx1 == 0

    block_idx1 = 0 if hc_left1 else np.max(np.where(np.array([b[0] for b in cv1.blocks]) < hc_idx1)[0])
    block1 = cv1.blocks[block_idx1][1]
    cv1.pos = int(block1[0]) if hc_left1 else int(block1[1])

    hc_idx2 = [idx for idx, clip in enumerate(r2.cigar) if clip[0] == CIGAR['hard-clip']][0]
    hc_left2 = hc_idx2 == 0

    block_idx2 = 0 if hc_left2 else np.max(np.where(np.array([b[0] for b in cv2.blocks]) < hc_idx2)[0])
    block2 = cv2.blocks[block_idx2][1]
    cv2.pos = int(block2[0]) if hc_left2 else int(block2[1])

    #TODO: handle inserted sequence btw fusion
    varseq1, varseq2 = '', ''
    refseq1 = r1.get_reference_sequence()
    refseq2 = r2.get_reference_sequence()
    bracket_dir1 = '[' if r1.is_reverse == r2.is_reverse else ']'
    bracket_dir2 = ']' if r1.is_reverse == r2.is_reverse else ']'
    cv1.ref = refseq1[0] if hc_left1 else refseq1[-1]
    cv2.ref = refseq2[0] if hc_left2 else refseq1[-1]
    if r1.is_reverse == r2.is_reverse:
        cv1.alt = '%s]%s:%d]%s' % (varseq1, cv2.chrom, cv2.pos-1, cv1.ref) if hc_left1 else \
                  '%s[%s:%d[%s' % (cv1.ref, cv2.chrom, cv2.pos+1, varseq1)
        cv2.alt = '%s]%s:%d]%s' % (varseq2, cv1.chrom, cv1.pos-1, cv2.ref) if hc_left2 else \
                  '%s[%s:%d[%s' % (cv2.ref, cv1.chrom, cv1.pos+1, varseq2)
    else:
        # contigs align on opposite strands
        cv1.alt = '%s[%s:%d[%s' % (varseq1, cv2.chrom, cv2.pos-1, cv1.ref) if hc_left1 else \
                  '%s]%s:%d]%s' % (cv1.ref, cv2.chrom, cv2.pos+1, varseq1)
        cv2.alt = '%s[%s:%d[%s' % (varseq2, cv1.chrom, cv1.pos-1, cv2.ref) if hc_left2 else \
                  '%s]%s:%d]%s' % (cv2.ref, cv1.chrom, cv1.pos+1, varseq2)

    cv1.vid = get_next_id(r1.query_name)
    cv2.vid = get_next_id(r2.query_name)
    cv1.parid, cv2.parid = cv2.vid, cv1.vid

    print(cv1.vcf_output())
    print(cv2.vcf_output())
    outbam.write(r1)
    outbam.write(r2)
    CrypticVariant.write_contig_info(ci_file, cv1, cv2)

    annotate_single_read(args, r1, juncs, ex_ref, ref_trees, genes=cv1.genes)
    annotate_single_read(args, r2, juncs, ex_ref, ref_trees, genes=cv2.genes)

def annotate_juncs(cv, read, locs, novel_juncs, ci_file):
    for junc in novel_juncs:
        pos1, pos2 = int(junc[1]), int(junc[2])
        junc_idx = [idx for idx, block in cv.blocks if block[1] == pos1][0]
        junc_type = read.cigar[junc_idx+1][0]
        if junc_type in GAPS or junc_type == CIGAR['soft-clip']:
            continue

        cp = CrypticVariant().from_read(read) # partner variant
        cp.genes = cv.genes
        varseq, refseq = '', read.get_reference_sequence()
        cpos = sum([v for c,v in read.cigar[:(junc_idx+1)] if c in AFFECT_CONTIG])
        rpos = sum([v for c,v in read.cigar[:(junc_idx+1)] if c in AFFECT_REF])
        cv.cpos, cp.cpos = cpos, cpos
        cv.pos, cp.pos = pos1-1, pos2+1

        cv.ref, cp.ref = refseq[rpos-1], refseq[rpos]
        cv.alt = '%s[%s:%d[%s' % (cv.ref, cv.chrom, cv.pos, varseq)
        cp.alt = '%s]%s:%d]%s' % (varseq, cp.chrom, cp.pos, cp.ref)

        cv.vid = get_next_id(read.query_name)
        cp.vid = get_next_id(read.query_name)
        cv.parid, cp.parid = cp.vid, cv.vid

        loc_left = '%s:%d' % (cv.chrom, pos1)
        loc_right = '%s:%d' % (cv.chrom, pos2)
        if not (loc_left in locs) and not (loc_right in locs):
            # neither end annotated, novel exon junction
            cv.cvtype, cp.cvtype = 'NEJ', 'NEJ'
        elif not (loc_left in locs and loc_right in locs):
            # one end unannotated, partial novel junction
            cv.cvtype, cp.cvtype = 'PNJ', 'PNJ'
        else:
            # both ends annotated, alternative splice site
            cv.cvtype, cp.cvtype = 'AS', 'AS'

        print(cv.vcf_output())
        print(cp.vcf_output())
        CrypticVariant.write_contig_info(ci_file, cv, cp)

def get_tx_juncs(read):
    tx_juncs = []
    unknown_juncs = []
    starts, ends = zip(*read.get_blocks())

    # merge adjacent 'junctions' (i.e. insertions)
    juncs = IntervalTree()
    for s, e in zip(starts, ends):
        juncs.addi(s, e)
    juncs.merge_overlaps(strict=False)
    starts = [junc[0] for junc in juncs]
    ends = [junc[1] for junc in juncs]

    chroms = [read.reference_name] * (len(starts)-1)
    tx_juncs = list(zip(chroms, ends[:-1], starts[1:]))
    tx_juncs = [junc for junc in tx_juncs if (junc[2] - junc[1]) > GAP_MIN]
    return tx_juncs

def annotate_single_read(args, read, juncs, ex_ref, ref_trees, outbam=None, genes=''):
    '''
    Annotate insertions, deletions and soft-clips on a single read
    '''
    ci_file = args.contig_info_file
    genes = get_overlapping_genes(read, ref_trees) if genes == '' else genes
    fusion = any([op == CIGAR['hard-clip'] and val >= CLIP_MIN for op, val in read.cigar])
    if genes == '' and not fusion:
        logging.info('No gene(s) intersecting read %s; skipping' % read.query_name)
        return

    # check for contig gaps or soft-clips
    has_gaps = any([op in GAPS and val >= GAP_MIN for op, val in read.cigar])
    has_scs = any([op == CIGAR['soft-clip'] and val >= CLIP_MIN for op, val in read.cigar])

    # check junctions
    tx_juncs = get_tx_juncs(read)
    unknown_juncs = ['%s:%s-%s' % (c, s, e) not in juncs[0] for c, s, e in tx_juncs]
    has_novel_juncs = any(unknown_juncs)

    # check for novel blocks
    chr_ref = get_chr_ref(read, ex_ref)
    has_novel_blocks = any([is_novel_block(block, chr_ref) for block in read.get_blocks()])

    if has_gaps or has_scs or has_novel_juncs or has_novel_blocks:
        cv = CrypticVariant().from_read(read)
        cv.genes = genes
        if has_gaps:
            annotate_gaps(cv, read, ci_file)
        if has_scs:
            annotate_softclips(cv, read, ci_file)
        if has_novel_juncs:
            novel_juncs = [list(x) for x in np.array(tx_juncs)[unknown_juncs]]
            annotate_juncs(cv, read, juncs[1], novel_juncs, ci_file)
        if has_novel_blocks:
            annotate_blocks(cv, read, chr_ref, ci_file)
        if outbam:
            outbam.write(read)
    else:
        logging.info('Nothing to annotate for read %s (read matches reference)' % read.query_name)

def annotate_contigs(args):
    '''
    Extract aligned contigs from supplied bam file and output
    annotated contig if it contains any novel bits
    '''
    ref_trees, ex_tree, ex_ref = get_gene_lookup(args.tx_ref_file)
    juncs = get_junc_lookup(args.junc_file)

    bam = pysam.AlignmentFile(args.bam_file, 'rc')
    bam_idx = pysam.IndexedReads(bam, multiple_iterators=True)
    bam_idx.build()
    outbam_file_unsort = '%s_unsorted.bam' % os.path.splitext(args.output_bam)[0]
    outbam = pysam.AlignmentFile(outbam_file_unsort, 'wb', template=bam)
    ci_file = args.contig_info_file

    logging.info('Checking contigs for non-reference content...')
    for read in bam.fetch(multiple_iterators=True):
        if read.reference_id < 0:
            logging.info('Skipping unmapped contig %s' % read.query_name)
            continue

        if read.query_name in record:
            # we have processed this read already (likely a read's partner)
            continue

        # only consider the contig if at least match_min bases align
        # to reference and at least match_perc_min of the read aligns
        rlen = read.reference_length
        qlen = float(read.query_length)
        if (rlen < MATCH_MIN) or (rlen / qlen) < MATCH_PERC_MIN:
            logging.info('Skipping contig %s: not enough bases match reference' % read.query_name)
            continue

        if all([op == CIGAR['match'] for op, val in read.cigar]):
            logging.info('Skipping contig %s: perfect match to reference' % read.query_name)
            continue

        is_hardclipped = any([op == CIGAR['hard-clip'] and val >= CLIP_MIN for op, val in read.cigar])
        if is_hardclipped:
            annotate_fusion(args, read, juncs, bam_idx, ex_ref, ref_trees, outbam)
        else:
            annotate_single_read(args, read, juncs, ex_ref, ref_trees, outbam)

    bam.close()
    outbam.close()

    # convert output sam file to bam, sort and index
    pysam.sort('-o', args.output_bam, outbam_file_unsort)
    pysam.index(args.output_bam)
    os.remove(outbam_file_unsort)

def main():
    args = parse_args()
    init_logging(args.log)
    print(VCF.get_header(args.sample))
    CrypticVariant.write_contig_header(args.contig_info_file)
    annotate_contigs(args)

if __name__ == '__main__':
    main()
