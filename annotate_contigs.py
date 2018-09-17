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
import re
import logging
import sys
from argparse import ArgumentParser
from intervaltree import Interval, IntervalTree
import ipdb

EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
PROGRAM_NAME = "annotate_contigs"

# CUT-OFF parameters
GAP_MIN = 7
CLIP_MIN = 30
MATCH_MIN = 30
MATCH_PERC_MIN = 0.3
INFO_STRING = "CID:ECN:CLEN:CSTRAND:CCIGAR:VSIZE:CSIZE:CVTYPE:GENES:PARID:ECC:AI:PVAL:CVQ"

# CIGAR specification codes
CIGAR = {'match': 0,
         'insertion': 1,
         'deletion': 2,
         'skipped': 3,
         'soft-clip': 4,
         'hard-clip': 5,
         'silent_deletion': 6}
GAPS = {1: 'insertion',
        2: 'deletion',
        6: 'silent_deletion'}
CLIPS = {4: 'soft',
         5: 'hard'}
GAPS_REF_ONLY = {2: 'deletion',
                 3: 'skipped',
                 5: 'hard-clip'} # criteria that signify gaps in the reference

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
    parser.add_argument(dest='sam_file',
                        metavar='SAM_FILE',
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

def get_gene_lookup(tx_ref_file):
    #TODO: Cache this step
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

def get_overlapping_gene(read):
    try:
        ref_tree = ref_trees[read.reference_name]
    except KeyError:
       logging.info('WARNING: reference chromosome %s (from read %s) not found in supplied reference' % (read.reference_name, read.query_name))
       return ''

    blocks = read.get_blocks()
    genes = []
    for block in blocks:
        bstart, bend = block
        gns = ref_tree.search(bstart, bend)
        genes.extend([gn[2] for gn in gns])
    genes = np.unique(np.array(genes))

    return('|'.join(genes))

class CrypticContig(object):
    '''
    Cryptic contig object
    '''
    def __init__(self, **kwargs):
        "Build an empty CrypticContig object"
        defaults = {
            "chrom": 'NA',
            "pos": -1,
            "strand": '.',
            "cid": '',
            "ref": '',
            "alt": '',
            "qual": -1,
            "cfilter": '.',
            "ecn": '',
            "clen": 0,
            "cpos": [-1,-1],
            "ccigar": '',
            "vsize": 0,
            "cvtype": '',
            "genes": '',
            "parid": '',
            "ecc": '',
            "ai": [-1, -1, -1],
            "pval": [-1, -1],
            "cvq": -1,
            "gap_idxs": [],
            "clip_idxs": [],
            "novel_juncs": [],
            "blocks": [],
            "cigar": [],
            "qseq": ''
        }
        for (prop, default) in defaults.items():
            setattr(self, prop, kwargs.get(prop, default))

    def from_read(self, read, novel_juncs):
        match_idxs = [idx for idx,cig in enumerate(read.cigar) if cig[0] == 0]
        self.chrom = read.reference_name
        self.strand = '-' if read.is_reverse else '+'
        self.blocks = [b for b in zip(match_idxs, read.get_blocks())]
        self.cigar = read.cigar
        self.cid = read.query_name
        self.ccigar = read.cigarstring
        self.novel_juncs = novel_juncs
        self.gap_idxs = [idx for idx, gap in enumerate(read.cigar) if gap[0] in GAPS and gap[1] >= GAP_MIN]
        self.clip_idxs = [idx for idx, clip in enumerate(read.cigar) if clip[0] in CLIPS and clip[1] >= GAP_MIN]
        self.clen =  sum([v for c,v in read.cigar if c in [0, 1, 4, 5]]) # count if match, insertion or clip TODO: parametrise
        self.qseq = read.query_sequence
        return self

    def get_info(self):
        cid = str(self.cid)
        ecn = str(self.ecn)
        clen = str(self.clen)
        cpos = str(self.cpos)
        cstrand = str(self.cstrand)
        ccigar = str(self.ccigar)
        vsize = str(self.vsize)
        cvtype = str(self.cvtype)
        genes = str(self.genes)
        parid = str(self.parid)
        ecc = str(self.ecc)
        ai = ','.join([str(a) for a in ai])
        pval = ','.join([str(pv) for pv in pval])
        cvq = str(self.cvq)
        return ':'.join([cid, ecn, clen, cpos, cstrand, ccigar, vsize, cvtype, genes, parid, ecc, ai, pval, cvq])

    def vcf_output(self):
        chrom = str(self.chrom)
        pos = str(self.pos)
        cid = str(self.cid)
        ref = str(self.ref)
        alt = str(self.alt)
        qual = str(self.qual)
        cfilter = str(self.cfilter)
        info_string = INFO_STRING
        info = self.get_info()
        return "\t".join([chrom, pos, cid, ref, alt, qual, cfilter, info_string, info])

    def annotate_gaps(self):
        '''
        Annotate gaps
        '''
        ipdb.set_trace()
        for gap_idx in self.gap_idxs:
            cigar = self.cigar[gap_idx]
            gtype, size = GAPS[cigar[0]], int(cigar[1])

            block_idx = 0 if gap_idx == 0 else np.max(np.where(np.array([b[0] for b in self.blocks]) < gap_idx)[0])
            block = self.blocks[block_idx][1]
            self.pos = int(block[1])

            # position of variant on contig
            cpos1 = sum([v for c,v in self.cigar[:gap_idx]])
            csize = size if self.cigar[gap_idx][0] == CIGAR['insertion'] else 0 # only insertions affect contig pos
            cpos2 = cpos1 + csize

            var_seq = ''
            if csize > 0:
                seq_pos1 = sum([v for c,v in self.cigar[:gap_idx] if c not in GAPS_REF_ONLY])
                seq_pos2 = seq_pos1 + csize
                var_seq = self.qseq[seq_pos1:seq_pos2]
            print(self.vcf_output())

    def annotate_clips(self):
        #TODO
        pass

    def annotate_juncs(self):
        #TODO
        pass

    def annotate_all(self, juncs, locs, ref_trees, ex_ref):
        if len(self.gap_idxs) > 0:
            self.annotate_gaps()
        if len(self.clip_idxs) > 0:
            self.annotate_clips()
        if len(self.novel_juncs) > 0:
            self.annotate_juncs()

def get_tx_juncs(read):
    tx_juncs = []
    unknown_juncs = []
    starts, ends = zip(*read.get_blocks())
    chroms = [read.reference_name] * (len(starts)-1)
    tx_juncs = list(zip(chroms, ends[:-1], starts[1:]))
    tx_juncs = [junc for junc in tx_juncs if (junc[2] - junc[1]) > GAP_MIN]
    return tx_juncs

def annotate_contigs(args):
    '''
    Extract aligned contigs from supplied sam file and output
    annotated contig if it contains any novel bits
    '''
    ref_trees, ex_ref = get_gene_lookup(args.tx_ref_file)
    juncs, locs = get_junc_lookup(args.junc_file)

    sam = pysam.AlignmentFile(args.sam_file, 'rc')
    outbam = pysam.AlignmentFile(args.output_bam, 'wb', template=sam)
    outbam_file_unsort = '%s_unsorted.bam' % os.path.splitext(args.output_bam)[0]

    logging.info('Checking contigs for non-reference content...')
    novel_contigs = []
    for read in sam.fetch():
        if read.reference_id < 0:
            logging.info('Skipping unmapped contig %s' % read.query_name)
            continue

        # only consider the contig if at least match_min bases align
        # to reference and at least match_perc_min of the read aligns
        rlen = read.reference_length
        qlen = float(read.query_length)
        if (rlen < MATCH_MIN) or (rlen / qlen) < MATCH_PERC_MIN:
            logging.info('Skipping contig %s: not enough bases match reference' % read.query_name)
            continue

        # check for contig gaps or clips
        has_gaps = any([op in GAPS and val >= GAP_MIN for op, val in read.cigar])
        has_clips = any([op in CLIPS and val >= CLIP_MIN for op, val in read.cigar])

        # check junctions
        tx_juncs = get_tx_juncs(read)
        unknown_juncs = ['%s:%s-%s' % (c, s, e) not in juncs for c, s, e in tx_juncs]
        has_novel_juncs = any(unknown_juncs)

        if has_gaps or has_clips or has_novel_juncs:
            novel_juncs = [list(x) for x in np.array(tx_juncs)[unknown_juncs]]
            cc = CrypticContig().from_read(read, novel_juncs)
            cc.annotate_all(juncs, locs, ref_trees, ex_ref)
            outbam.write(read)

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
