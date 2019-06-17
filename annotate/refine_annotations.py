'''
Module      : refine_annotations
Description : Perform further filtering on annotated contigs
Copyright   : (c) Marek Cmero, Feb 2019
License     : TBD
Maintainer  : MAREK.CMERO@MCRI.EDU.AU
Portability : POSIX
'''

import numpy as np
import pandas as pd
import re
import sys
import logging
import pysam
import bedtool_helper as bh
import make_supertranscript as ms
import annotate_contigs as ac
import pybedtools as pbt
import constants
from pybedtools import BedTool
from argparse import ArgumentParser
from utils import init_logging, exit_with_error

PROGRAM_NAME = 'refine_annotations'

SPLICE_VARS = ['AS']
SV_VARS = ['DEL', 'INS', 'UN']
NOVEL_BLOCKS = ['EE', 'NE']
NOVEL_JUNCS = ['PNJ', 'NEJ']
FUSIONS = ['FUS']

def parse_args():
    '''
    Parse command line arguments.
    Returns Options object with command line argument values as attributes.
    Will exit the program on a command line error.
    '''
    description = 'Filter contigs'
    parser = ArgumentParser(description=description)
    parser.add_argument('--log',
                        metavar='LOG_FILE',
                        type=str,
                        help='record program progress in LOG_FILE')
    parser.add_argument(dest='contig_info_file',
                        metavar='CONTIG_INFO_FILE',
                        type=str,
                        help='''Contig info file''')
    parser.add_argument(dest='vcf_file',
                        metavar='VCF_FILE',
                        type=str,
                        help='''Contig VCF file''')
    parser.add_argument(dest='bam_file',
                        metavar='BAM_FILE',
                        type=str,
                        help='''Contig BAM file''')
    parser.add_argument(dest='tx_ref_file',
                        type=str,
                        metavar='TX_REF_FILE',
                        help='''Transcriptiome GTF reference file.''')
    parser.add_argument(dest='fasta',
                        metavar='FASTA',
                        type=str,
                        help='''Genome fasta file''')
    parser.add_argument(dest='contig_out_file',
                        metavar='CONTIG_OUT_FILE',
                        type=str,
                        help='''Contig tsv output file''')
    parser.add_argument(dest='bam_out_file',
                        metavar='BAM_OUT_FILE',
                        type=str,
                        help='''Contig BAM output file''')
    parser.add_argument('--minClip',
                        metavar='MIN_CLIP',
                        type=int,
                        help='''Minimum novel block or softclip size.''')
    parser.add_argument('--minGap',
                        metavar='MIN_GAP',
                        type=int,
                        help='''Minimum gap (deletion or insertion) size.''')
    return parser.parse_args()

def set_globals(args):
    global MIN_CLIP
    global MIN_GAP

    if args.minClip:
        MIN_CLIP = args.minClip
    else:
        MIN_CLIP = constants.DEFAULT_MIN_CLIP

    if args.minGap:
        MIN_GAP = args.minGap
    else:
        MIN_GAP = constants.DEFAULT_MIN_GAP

def is_valid_motif(left_id, right_id, block_seqs):
    try:
        if left_id != '' and right_id != '':
            valid_sense = block_seqs[left_id] == 'AG' and block_seqs[right_id] == 'GT'
            valid_antisense = block_seqs[left_id] == 'AC' and block_seqs[right_id] == 'CT'
            return valid_sense or valid_antisense
        elif left_id != '':
            return block_seqs[left_id] in ['AG', 'AC']
        elif right_id != '':
            return block_seqs[right_id] in ['GT', 'CT']
        else:
            return False
    except KeyError:
        # occurs if chrom not in reference
        logging.info('''WARNING: one of the following motif Locations could not be
                        retrieved from the provided reference: %s, %s''' % (left_id, right_id))
        return False

def get_valid_motif_vars(variants, args):
    '''
    Checks whether the given variant has a canonical splice site motif,
    returning those variant IDs with valid motifs. Novel exons are checked
    for valid donor and acceptor sites (in either transcriptional direction),
    while extended exons on the left and right are checked for either valid
    donor or acceptor sites (in either transcriptional direction).
    '''
    # get VCF data of given variants
    vcf = ms.load_vcf_file(args.vcf_file)
    vcf = vcf[vcf[2].apply(lambda x: x in variants.variant_id.values)]
    vcf = vcf[vcf[3].apply(lambda x: np.invert(pd.isnull(x)))]

    # construct motif locations for which to extract sequence
    b_blocks = vcf[vcf[4].apply(lambda x: x.startswith('[') and x.endswith(']'))]
    r_blocks = vcf[vcf[4].apply(lambda x: x.startswith(']'))]
    l_blocks = vcf[vcf[4].apply(lambda x: x.endswith('['))]
    left = pd.concat([l_blocks, b_blocks])
    right = pd.concat([r_blocks, b_blocks])
    mlocs = pd.DataFrame({'chr': pd.concat([left[0], right[0]]),
                          'start': pd.concat([left[1] - 3, right[1] + right[3].apply(len)-1])})
    mlocs['end'] = mlocs['start'] + 2
    mlocs = mlocs.drop_duplicates()

    # ensure coords aren't negative or exceed chrom len
    mlocs = mlocs[mlocs['start'] >= 0]
    chr_sizes = pbt.chromsizes('hg38') #TODO: allow reference to be argument
    for chrom in np.unique(mlocs['chr'].values):
        ref_chrom = 'chr%s' % chrom if chrom != 'MT' else 'chrM'
        try:
            chr_max = chr_sizes[ref_chrom][1]
            over_limit = np.logical_or(mlocs.start > chr_max, mlocs.end > chr_max)
            mlocs = mlocs.drop(mlocs[np.logical_and(mlocs['chr'] == chrom, over_limit)].index)
        except KeyError:
            logging.info("WARNING: chrom %s does not exist in hg38 reference." % ref_chrom)
            continue

    # extract sequences
    g = BedTool.from_dataframe(mlocs).remove_invalid()
    g = g.sequence(fi=args.fasta)
    bs = bh.get_block_seqs(g)

    # check whether variants have valid motifs
    # left blocks
    valid_vars = []
    left_ids = ['%s:%d-%d' % loc for loc in zip(l_blocks[0], l_blocks[1]-3, l_blocks[1]-1)]
    valid_left = [is_valid_motif(lid, '', bs) for lid in left_ids]
    valid_vars.extend(l_blocks[valid_left][2].values)

    # right blocks
    rpos = r_blocks[1] + r_blocks[3].apply(len)-1
    right_ids = ['%s:%d-%d' % loc for loc in zip(r_blocks[0], rpos, rpos+2)]
    valid_right = [is_valid_motif('', rid, bs) for rid in right_ids]
    valid_vars.extend(r_blocks[valid_right][2].values)

    # both blocks (novel exons)
    rpos = b_blocks[1] + b_blocks[3].apply(len)-1
    left_ids = ['%s:%d-%d' % loc for loc in zip(b_blocks[0], b_blocks[1]-3, b_blocks[1]-1)]
    right_ids = ['%s:%d-%d' % loc for loc in zip(b_blocks[0], rpos, rpos+2)]
    valid_both = [is_valid_motif(lid, rid, bs) for lid, rid in zip(left_ids, right_ids)]
    valid_vars.extend(b_blocks[valid_both][2].values)

    return np.unique(valid_vars)

def check_overlap(ex_trees, chrom, start, end, check_size):
    '''
    Checks whether variant overlaps an exonic region.
    For deletions, at least MIN_GAP bp of the deletion
    must be within the exon body.
    '''
    olap = False
    try:
        olap = ex_trees[chrom].overlaps(start, end)
        if olap and check_size:
            olap_se = ex_trees[chrom].overlap(start, end)
            es, ee = [(x[0], x[1]) for x in olap_se][0]
            size = min([ee, end]) - start if start >= es \
                                          else end - max([es, start])
            olap = size >= MIN_GAP
    except KeyError:
        logging.info('WARNING: chrom %s not found in provided reference.' % chrom)
    return olap

def overlaps_exon(sv, ex_trees):
    pos1 = sv['pos1'].split(':')
    pos2 = sv['pos2'].split(':')

    chrom = pos1[0]
    start = int(pos1[1].split('(')[0])
    end = int(pos2[1].split('(')[0])
    end = start + 1 if sv['variant_type'] != 'DEL' else end

    check_size = sv['variant_type'] == 'DEL'
    olap = check_overlap(ex_trees, chrom, start, end, check_size)
    if sv['variant_type'] == 'FUS':
        chrom = pos2[0]
        start = int(pos2[1].split('(')[0])
        olap = olap or check_overlap(ex_trees, chrom, start,
                                     start + 1, check_size)

    return olap

def get_contigs_to_keep(args):
    '''
    Return contigs matching criteria:
    - novel block size > MIN_CLIP
    - novel blocks are spliced in some way
    - novel exons have corresponding novel splice
      sites and novel splice donor/acceptor motifs
    - TSVs occur in exonic regions
    '''
    try:
        cinfo_file = args.contig_info_file
        contigs = pd.read_csv(cinfo_file, sep='\t')
    except IOError as exception:
        exit_with_error(str(exception), constants.EXIT_FILE_IO_ERROR)

    # general tests - var size and splicing
    contigs['large_varsize'] = contigs.contig_varsize > MIN_CLIP
    contigs['is_contig_spliced'] = contigs.contig_cigar.apply(lambda x: bool(re.search('N', x)))

    # test for valid splice motifs
    contigs['valid_motif'] = False
    check_motifs = contigs.variant_type.isin(NOVEL_JUNCS \
                                             + NOVEL_BLOCKS \
                                             + ['DEL'])
    if any(check_motifs.values):
        valid_motif_vars = get_valid_motif_vars(contigs[check_motifs], args)
        contigs.loc[check_motifs, 'valid_motif'] = \
            contigs.variant_id.isin(valid_motif_vars)

    # check whether exons contain matching novel/partial novel junctions
    spliced_exons = []
    exons = contigs[contigs.variant_type.isin(NOVEL_BLOCKS)]
    novel_juncs = contigs[contigs.variant_type.isin(NOVEL_JUNCS)]
    for idx,row in exons.iterrows():
        back_junc = novel_juncs.pos2 == row['pos1']
        front_junc = novel_juncs.pos1 == row['pos2']
        matching_juncs = novel_juncs[np.logical_or(back_junc, front_junc)]
        if len(matching_juncs) > 0:
            spliced_exons.append(row['variant_id'])
    contigs['spliced_exon'] = contigs.variant_id.isin(spliced_exons)

    # novel exon contigs (spliced, valid motif and large variant size)
    is_novel_exon = np.logical_and(contigs.spliced_exon, contigs.valid_motif)
    is_novel_exon = np.logical_and(is_novel_exon, contigs.large_varsize)
    ne_vars = contigs[is_novel_exon].variant_id.values

    # check whether TSVs meets size requirements
    is_sv = contigs.variant_type.isin(SV_VARS)
    is_fus = contigs.variant_type.isin(FUSIONS)
    large_clip = np.logical_or(contigs.varsize >= MIN_GAP,
                               contigs.contig_varsize >= MIN_GAP)
    keep_sv = np.logical_and(large_clip, is_sv)
    keep_sv = np.logical_or(keep_sv, is_fus)

    # check whether TSVs are within exons
    contigs['overlaps_exon'] = False
    exonic_var = contigs.variant_type.isin(['PNJ', 'EE', 'AS'])
    contigs.loc[exonic_var, 'overlaps_exon'] = True
    gene_tree, ex_trees, ex_ref = ac.get_gene_lookup(args.tx_ref_file)
    contigs.loc[keep_sv, 'overlaps_exon'] = contigs[keep_sv].apply(overlaps_exon,
                                                                   axis=1, args=(ex_trees,))
    sv_vars = contigs[np.logical_and(keep_sv, contigs.overlaps_exon)].variant_id.values

    # keep all splice vars
    as_vars = contigs.variant_id.values[contigs.variant_type.isin(SPLICE_VARS)]

    # ensure retained introns are spliced in some way (to distinguish from pre-mRNAs)
    retained_intron = contigs.variant_type == 'RI'
    ri_vars = contigs[np.logical_and(retained_intron, contigs.large_varsize)].variant_id.values

    # collate contigs to keep
    keep_vars = np.unique(np.concatenate([ri_vars, as_vars, ne_vars, sv_vars]))
    contigs['variant_of_interest'] = contigs.variant_id.isin(keep_vars)
    keep_contigs = contigs[contigs.variant_of_interest].contig_id.values
    contigs = contigs[contigs.contig_id.isin(keep_contigs)]
    contigs.to_csv(args.contig_out_file, sep='\t', index=None)

    return(keep_contigs)

def write_output(args, keep_contigs):
    cvars_file = args.vcf_file
    try:
        vcf = pd.read_csv(cvars_file, sep='\t', header=None, comment='#')
        cvf = open(cvars_file, 'r')
        for line in cvf:
            if not line.startswith('#'):
                break
            print(line.strip())
        cvf.close()
    except IOError as exception:
        exit_with_error(str(exception), constants.EXIT_FILE_IO_ERROR)

    vcf = vcf[vcf[7].apply(lambda x: x.split(';')[0].split('=')[1] in keep_contigs)]
    vcf.to_csv(sys.stdout, sep='\t', index=False, header=False)

def write_bam(args, keep_contigs):
    bam_file = args.bam_file
    bam = pysam.AlignmentFile(bam_file, 'rb')
    outbam = pysam.AlignmentFile(args.bam_out_file, 'wb', template=bam)

    for read in bam.fetch():
        if read.query_name in keep_contigs:
           outbam.write(read)

def main():
    args = parse_args()
    init_logging(args.log)
    set_globals(args)
    keep_contigs = get_contigs_to_keep(args)
    write_output(args, keep_contigs)
    write_bam(args, keep_contigs)

if __name__ == '__main__':
    main()
