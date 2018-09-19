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
from cv_vcf import CrypticVariant
import ipdb

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
GAPS = [CIGAR['insertion'], CIGAR['deletion'], CIGAR['silent_deletion']]
# any cigar criteria that is >0 bp on an aligned contig
AFFECT_CONTIG = [CIGAR['insertion'], CIGAR['match'], CIGAR['soft-clip'], CIGAR['hard-clip']]

def cached(cachefile):
    '''
    source: https://datascience.blog.wzb.eu/2016/08/12/a-tip-for-the-impatient-simple-caching-with-python-pickle-and-decorators/
    A function that creates a decorator which will use "cachefile"
    for caching the results of the decorated function "fn".
    '''
    def decorator(fn):  # define a decorator for a function "fn"
        def wrapped(*args, **kwargs):   # define a wrapper that will finally call "fn" with all arguments
          # if cache exists -> load it and return its content
          if os.path.exists(cachefile):
              with open(cachefile, 'rb') as cachehandle:
                print("using cached result from '%s'" % cachefile)
                return pickle.load(cachehandle)

          # execute the function with all arguments passed
          res = fn(*args, **kwargs)

          # write to cache file
          with open(cachefile, 'wb') as cachehandle:
            print("saving result to cache '%s'" % cachefile)
            pickle.dump(res, cachehandle)

          return res

        return wrapped

    return decorator   # return this "customized" decorator that uses "cachefile"

def exit_with_error(message, exit_status):
    '''
    Print an error message to stderr, prefixed by the program name and 'ERROR'.
    Then exit program with supplied exit status.
    Arguments:
        message: an error message as a string.
        exit_status: a positive integer representing the exit status of the
            program.
    '''
    logging.error(message)
    print("{} ERROR: {}, exiting".format(PROGRAM_NAME, message), file=sys.stderr)
    sys.exit(exit_status)

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
    parser.add_argument(dest='bam_file',
                        metavar='BAM_FILE',
                        type=str,
                        help='''SAM or BAM format file containing contig alignments''')
    parser.add_argument(dest='output_bam',
                        metavar='OUTPUT_BAM',
                        type=str,
                        help='''BAM file to write contigs which pass filtering''')
    parser.add_argument(dest='junc_file',
                        metavar='JUNC_FILE',
                        type=str,
                        help='''Reference file containing transcripts and their
                        respective splice junctions.''')
    parser.add_argument(dest='tx_ref_file',
                        type=str,
                        metavar='TX_REF_FILE',
                        help='''Transcriptiome GTF reference file.''')
    parser.add_argument(dest='vcf_header_file',
                        type=str,
                        metavar='VCF_HEADER_FILE',
                        help='''VCF header file for cryptic variant annotation.''')
    return parser.parse_args()

def init_logging(log_filename):
    '''
    If the log_filename is defined, then
    initialise the logging facility, and write log statement
    indicating the program has started, and also write out the
    command line from sys.argv
    Arguments:
        log_filename: either None, if logging is not required, or the
            string name of the log file to write to
    Result:
        None
    '''
    if log_filename is not None:
        logging.basicConfig(filename=log_filename,
                            level=logging.DEBUG,
                            filemode='w',
                            format='%(asctime)s %(levelname)s - %(message)s',
                            datefmt='%m-%d-%Y %H:%M:%S')
        logging.info('program started')
    logging.info('command line: %s', ' '.join(sys.argv))

@cached('gene_lookup_cache.pickle')
def get_gene_lookup(tx_ref_file):
    '''
    Generate start/end coordinate reference
    for genes and output as an interval tree
    dictionary. Also output dataframe containing
    chromosome, start and ends for all exons.
    '''
    ref_trees = None
    if tx_ref_file != '':
        logging.info('Generating lookup for genes...')
        tx_ref = pd.read_csv(tx_ref_file, comment='#', sep='\t', header=None)
        gn_ref = tx_ref[tx_ref[2] == 'gene']
        gn_ref.loc[:, 'gene'] = gn_ref[8].apply(lambda x: re.search('gene_name "([\w\-\.\/]+)"', x).group(1))
        gn_ref = gn_ref[[0, 3, 4, 'gene']]
        gn_ref.columns = ['chrom', 'start', 'end', 'gene']
        gn_ref = gn_ref.drop_duplicates()

        ref_trees = {}
        chroms = np.unique(gn_ref.chrom.values)
        for chrom in chroms:
            chr_ref = gn_ref[gn_ref.chrom == chrom]
            ref_tree = IntervalTree()
            for idx, row in chr_ref.iterrows():
                ref_tree[row.start:row.end] = row.gene
            ref_trees[chrom] = ref_tree

        ex_ref = tx_ref[tx_ref[2] == 'exon']
        ex_ref = ex_ref[[0, 3, 4]].drop_duplicates()
        ex_ref.columns = ['chrom', 'start', 'end']
        ex_ref.start = ex_ref.start - 1
    return ref_trees, ex_ref

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

    return juncs, locs

def get_overlapping_genes(read, ref_trees):
    try:
        ref_tree = ref_trees[read.reference_name]
    except KeyError:
       logging.info('WARNING: reference chromosome %s (from read %s) not found in supplied reference' %
                    (read.reference_name, read.query_name))
       return ''

    blocks = read.get_blocks()
    genes = []
    for block in blocks:
        bstart, bend = block
        gns = ref_tree.search(bstart, bend)
        genes.extend([gn[2] for gn in gns])
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

def get_next_id(qname, record):
    if qname not in record:
        vid = qname + 'a'
        record[qname] = [vid]
        return vid, record
    else:
        vid = qname + get_next_letter(record[qname][-1][-1])
        record[qname].append(vid)
        return vid, record

def annotate_gaps(cv, read, record):
    '''
    Annotate deletions and insertions
    '''
    gap_idxs = [idx for idx, gap in enumerate(read.cigar) if gap[0] in GAPS and gap[1] >= GAP_MIN]
    for gap_idx in gap_idxs:
        cigar = read.cigar[gap_idx]
        gtype, cv.vsize = GAPS[cigar[0]], int(cigar[1])

        block_idx = 0 if gap_idx == 0 else np.max(np.where(np.array([b[0] for b in cv.blocks]) < gap_idx)[0])
        block = cv.blocks[block_idx][1]
        cv.pos = int(block[1])

        # position of variant on contig
        cv.cpos[0] = sum([v for c,v in read.cigar[:gap_idx]])
        cv.cvsize = cv.vsize if read.cigar[gap_idx][0] == CIGAR['insertion'] else 0 # only insertions affect contig pos
        cv.cpos[1] = cv.cpos[0] + cv.cvsize

        var_seq = ''
        if cv.cvsize > 0:
            seq_pos1 = sum([v for c,v in read.cigar[:gap_idx] if c in AFFECT_CONTIG])
            seq_pos2 = seq_pos1 + cv.cvsize
            seq = read.query_sequence[(seq_pos1-1):seq_pos2]
            cv.ref, cv.alt = seq[:1], seq
        else:
            start, end = cv.cpos[0]-1, (cv.cpos[0] + cv.vsize)
            seq = read.get_reference_sequence()[start:end]
            cv.ref, cv.alt  = seq, seq[:1]

        cv.cvtype = 'DEL' if gtype in ['deletion', 'silent_deletion'] else 'INS'
        cv.vid, record = get_next_id(read.query_name, record)
        print(cv.vcf_output())
    return record

def annotate_softclips(cv, read, record):
    sc_idxs = [idx for idx, clip in enumerate(read.cigar) if clip[0] == CIGAR['soft-clip'] and clip[1] >= CLIP_MIN]
    for sc_idx in sc_idxs:
        cigar = read.cigar[sc_idx]
        cv.cvtype, cv.vsize, cv.cvsize = 'UT', int(cigar[1]), int(cigar[1])
        sc_left = sc_idx == 0

        block_idx = 0 if sc_idx == 0 else np.max(np.where(np.array([b[0] for b in cv.blocks]) < sc_idx)[0])
        block = cv.blocks[block_idx][1]
        cv.pos = int(block[0]) if sc_left else int(block[1])

        rcigar = read.cigar[::-1] if cv.strand == '-' else read.cigar
        cv.cpos[0] = sum([v for c,v in rcigar[:sc_idx] if c in AFFECT_CONTIG])
        cv.cpos[1] = cv.cpos[0] + cv.cvsize

        varseq = read.query_sequence[cv.cpos[0]:cv.cpos[1]]
        refseq = read.get_reference_sequence()
        cv.ref = refseq[0] if sc_left else refseq[-1]
        cv.alt = '%s]%s:%d]%s' % (varseq, cv.chrom, cv.pos-1, cv.ref) if sc_left else \
                 '%s[%s:%d[%s' % (cv.ref, cv.chrom, cv.pos+1, varseq)

        cv.vid, record = get_next_id(read.query_name, record)
        print(cv.vcf_output())
    return record

def annotate_block(cv, read, block, block_idx, ex_chr):
    '''
    For fusions with non-gene sequence at break
    '''
    block_size = block[1] - block[0]
    olapping = ex_chr[np.logical_and(block[0] < ex_chr.start, block[1] > ex_chr.end)]
    left = ex_chr[np.logical_and(block[1] > ex_chr.start, block[1] <= ex_chr.end)]
    right = ex_chr[np.logical_and(block[0] >= ex_chr.start, block[0] < ex_chr.end)]

    seq_pos1 = sum([e - s for s, e in read.get_blocks()[:block_idx]])
    seq_pos2 = seq_pos1 + block_size
    block_seq = read.query_sequence[seq_pos1:seq_pos2]
    varseq = ''

    if len(olapping) == 0 and len(right) == 0 and len(left) == 0:
        # block is entirely intronic or intergenic
        varseq = block_seq
    elif len(left) > 0 and len(right) > 0:
        # intron retention; split into two
        pass
    elif len(left) > 0 or len(right) > 0:
        # get variant sequence not matching reference
        olap_size = block[1] - min(left.start.values) if len(left)>0 else max(right.end.values) - block[0]
        varseq = block_seq[:-olap_size] if len(left)>0 else block_seq[olap_size:]
    elif len(olapping) > 0:
        # if block overlaps an exon, must split into two
        olap_left = min(olapping.start.values) - block[0]
        varseq = block_seq[:olap_left]
        vpos1 = size if clip_left else 0
        if read.is_reverse:
            vpos1 = contig_size - len(varseq)
            vpos1 = vpos1 if not clip_left else vpos1 - size
        vpos2 = vpos1 + len(varseq)
    return varseq

def annotate_fusion(read, juncs, locs, bam_idx, ex_ref, ref_trees, outbam, record):
    try:
        r1, r2 = bam_idx.find(read.query_name)
    except ValueError:
        logging.info('WARNING: found >2 reads matching hard-clipped read %s; cannot process' % read.query_name)
        return record

    cv1 = CrypticVariant().from_read(r1)
    cv2 = CrypticVariant().from_read(r2)
    cv1.cvtype, cv2.cvtype = 'FUS', 'FUS'

    cv1.genes = get_overlapping_genes(r1, ref_trees)
    cv2.genes = get_overlapping_genes(r2, ref_trees)
    if cv1.genes == '' and cv2.genes == '':
        # no intersecting gene, this is not an interesting fusion
        logging.info('No gene(s) intersecting candidate fusion contig %s; skipping' % read.query_name)
        record[read.query_name] = []
        return record

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
    cv1.ref = refseq1[0] if hc_left1 else refseq1[-1]
    cv1.alt = '%s]%s:%d]%s' % (varseq1, cv2.chrom, cv2.pos-1, cv1.ref) if hc_left1 else \
              '%s[%s:%d[%s' % (cv1.ref, cv2.chrom, cv2.pos+1, varseq1)
    cv2.ref = refseq2[0] if hc_left2 else refseq1[-1]
    cv2.alt = '%s]%s:%d]%s' % (varseq2, cv1.chrom, cv1.pos-1, cv2.ref) if hc_left2 else \
              '%s[%s:%d[%s' % (cv2.ref, cv1.chrom, cv1.pos+1, varseq2)

    cv1.vid, record = get_next_id(r1.query_name, record)
    cv2.vid, record = get_next_id(r2.query_name, record)
    cv1.parid, cv2.parid = cv2.vid, cv1.vid

    print(cv1.vcf_record())
    print(cv2.vcf_record())
    outbam.write(r1)
    outbam.write(r2)

    record = annotate_single_read(r1, juncs, locs, ex_ref, record)
    record = annotate_single_read(r2, juncs, locs, ex_ref, record)

    return record

def annotate_juncs(cv, read, locs, novel_juncs, record):
    for junc in novel_juncs:
        pos1, pos2 = int(junc[1]), int(junc[2])
        junc_idx = [idx for idx, block in cv.blocks if block[1] == pos1][0]
        junc_type = read.cigar[junc_idx+1][0]
        if junc_type in GAPS or junc_type == CIGAR['soft-clip']:
            continue

        cp = CrypticVariant().from_read(read) # partner variant
        cp.genes = cv.genes
        varseq, refseq = '', read.get_reference_sequence()
        cpos = sum([v for c,v in read.cigar[:(junc_idx+1)]])
        cv.cpos[0], cv.cpos[1] = cpos, cpos
        cp.cpos[0], cp.cpos[1] = cv.cpos[0], cv.cpos[1]
        cv.pos, cp.pos = pos1-1, pos2+1

        cv.ref, cp.ref = refseq[cpos-1], refseq[cpos+1]
        cv.alt = '%s[%s:%d[%s' % (cv.ref, cv.chrom, cv.pos, varseq)
        cp.alt = '%s]%s:%d]%s' % (varseq, cp.chrom, cp.pos, cp.ref)

        cv.vid, record = get_next_id(read.query_name, record)
        cp.vid, record = get_next_id(read.query_name, record)
        cv.parid, cp.parid = cp.vid, cv.vid

        # TODO: annotate variant type
        # TODO: test
        ipdb.set_trace()
        print(cv.vcf_output())
        print(cp.vcf_output())
    return record

def get_tx_juncs(read):
    tx_juncs = []
    unknown_juncs = []
    starts, ends = zip(*read.get_blocks())
    chroms = [read.reference_name] * (len(starts)-1)
    tx_juncs = list(zip(chroms, ends[:-1], starts[1:]))
    tx_juncs = [junc for junc in tx_juncs if (junc[2] - junc[1]) > GAP_MIN]
    return tx_juncs

def annotate_single_read(read, juncs, locs, ex_ref, record, outbam=None, ref_trees=None):
    '''
    Annotate insertions, deletions and soft-clips on a single read
    '''
    #TODO: add support for single block read retained intron
    genes = ''
    if ref_trees:
        genes = get_overlapping_genes(read, ref_trees)
        if genes == '':
            logging.info('No gene(s) intersecting read %s; skipping' % read.query_name)
            return record

    # check for contig gaps or soft-clips
    has_gaps = any([op in GAPS and val >= GAP_MIN for op, val in read.cigar])
    has_scs = any([op == CIGAR['soft-clip'] and val >= CLIP_MIN for op, val in read.cigar])

    # check junctions
    tx_juncs = get_tx_juncs(read)
    unknown_juncs = ['%s:%s-%s' % (c, s, e) not in juncs for c, s, e in tx_juncs]
    has_novel_juncs = any(unknown_juncs)

    if has_gaps or has_scs or has_novel_juncs:
        cv = CrypticVariant().from_read(read)
        cv.genes = genes
        if has_gaps:
            record = annotate_gaps(cv, read, record)
        if has_scs:
            record = annotate_softclips(cv, read, record)
        if has_novel_juncs:
            novel_juncs = [list(x) for x in np.array(tx_juncs)[unknown_juncs]]
            record = annotate_juncs(cv, read, locs, novel_juncs, record)
        #TODO: intron retention
        if outbam:
            outbam.write(read)
    else:
        logging.info('Nothing to annotate for read %s (read matches reference)' % read.query_name)

    return record

def annotate_contigs(args):
    '''
    Extract aligned contigs from supplied bam file and output
    annotated contig if it contains any novel bits
    '''
    ref_trees, ex_ref = get_gene_lookup(args.tx_ref_file)
    juncs, locs = get_junc_lookup(args.junc_file)

    bam = pysam.AlignmentFile(args.bam_file, 'rc')
    bam_idx = pysam.IndexedReads(bam)
    bam_idx.build()
    outbam = pysam.AlignmentFile(args.output_bam, 'wb', template=bam)
    outbam_file_unsort = '%s_unsorted.bam' % os.path.splitext(args.output_bam)[0]

    logging.info('Checking contigs for non-reference content...')
    novel_contigs = []
    record = {}
    for read in bam.fetch():
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

        is_hardclipped = any([op == CIGAR['hard-clip'] and val >= CLIP_MIN for op, val in read.cigar])
        if is_hardclipped:
            record = annotate_fusion(read, juncs, locs, bam_idx, ex_ref, ref_trees, outbam, record)
        else:
            record = annotate_single_read(read, juncs, locs, ex_ref, record, outbam, ref_trees)

    sam.close()
    outbam.close()

    # convert output sam file to bam, sort and index
    pysam.sort('-o', outbam, outbam_file_unsort)
    pysam.index(args.outbam_file)
    os.remove(outbam_file_unsort)

def main():
    args = parse_args()
    init_logging(args.log)
    with open(args.vcf_header_file, 'rb') as headerf:
        headerf.readlines()
    annotate_contigs(args)

if __name__ == '__main__':
    main()
