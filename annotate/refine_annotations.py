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
import ipdb
import bedtool_helper as bh
import make_supertranscript as ms
import annotate_contigs as ac
from pybedtools import BedTool
from argparse import ArgumentParser
from utils import init_logging, exit_with_error

EXIT_FILE_IO_ERROR = 1

MIN_NOVEL_EXON_SIZE = 20
SPLICE_VARS = ['AS']
SV_VARS = ['FUS', 'DEL', 'INS', 'UN']
NOVEL_BLOCKS = ['EE', 'NE']

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

    return parser.parse_args()

def is_valid_splice_motif(row, fasta):
    '''
    Checks whether the given variant has a canonical splice site motif.
    Novel exons are checked for valid donor and acceptor sites (in either
    transcriptional direction), while extended exons on the left and right
    are checked for either valid donor or acceptor sites (in either trans-
    criptional direction).
    '''
    both = row[4].startswith('[') and row[4].endswith(']')
    left = row[4].endswith('[') or both
    right = row[4].startswith(']') or both

    mloc = pd.DataFrame()
    slen = len(row[3])
    if left:
        df = {'chr': row[0], 'pos1': row[1] - 3, 'pos2': row[1] - 1}
        mloc = mloc.append(df, ignore_index=True)
    if right:
        df = {'chr': row[0], 'pos1': row[1] + slen, 'pos2': row[1] + slen + 2}
        mloc = mloc.append(df, ignore_index=True)
    mloc['pos1'], mloc['pos2'] = mloc.pos1.map(int), mloc.pos2.map(int)

    g = BedTool.from_dataframe(mloc)
    g = g.sequence(fi=fasta)
    bs = bh.get_block_seqs(g)
    if len(bs) == 0:
        return False

    left_loc = '%s:%d-%d' % (row[0], row[1] - 3, row[1] - 1)
    right_loc = '%s:%d-%d' % (row[0], row[1] + slen, row[1] + slen + 2)

    if both:
        valid_sense = bs[left_loc] == 'AG' and bs[right_loc] == 'GT'
        valid_antisense = bs[left_loc] == 'AC' and bs[right_loc] == 'CT'
        return valid_sense or valid_antisense
    elif left:
        return bs[left_loc] in ['AG', 'AC']
    elif right:
        return bs[right_loc] in ['GT', 'CT']

def get_valid_motif_vars(variants, args):
    '''
    Ensure novel exons contain valid canonical splice motifs
    '''
    fasta = args.fasta
    vcf = ms.load_vcf_file(args.vcf_file)
    vcf = vcf[vcf[2].apply(lambda x: x in variants.variant_id.values)]

    valid_motifs = vcf.apply(is_valid_splice_motif, axis=1, args=(fasta,))
    valid_vars = vcf[valid_motifs][2].values

    return valid_vars

def overlaps_exon(sv, ex_trees):
    pos1 = sv['pos1'].split(':')
    pos2 = sv['pos2'].split(':')

    chrom = pos1[0]
    start = int(pos1[1].split('(')[0])
    end = int(pos2[1].split('(')[0])
    end = start + 1 if sv['variant_type'] != 'DEL' else end

    olap = ex_trees[chrom].overlaps(start, end)
    if sv['variant_type'] == 'FUS':
        chrom = pos2[0]
        start = int(pos2[1].split('(')[0])
        olap = olap or ex_trees[chrom].overlaps(start, start+1)

    return olap

def get_contigs_to_keep(args):
    '''
    Return contigs matching criteria:
    - novel block size > MIN_NOVEL_EXON_SIZE
    - novel blocks are spliced in some way
    - novel exons have corresponding novel splice
      sites and novel splice donor/acceptor motifs
    - TSVs occur in exonic regions
    '''
    try:
        cinfo_file = args.contig_info_file
        contigs = pd.read_csv(cinfo_file, sep='\t')
    except IOError as exception:
        exit_with_error(str(exception), EXIT_FILE_IO_ERROR)

    # general tests - var size and splicing
    contigs['large_varsize'] = contigs.contig_varsize > MIN_NOVEL_EXON_SIZE
    contigs['is_contig_spliced'] = contigs.contig_cigar.apply(lambda x: bool(re.search('N', x)))

    # test for valid splice motifs
    contigs['valid_motif'] = False
    check_motifs = contigs.variant_type.apply(lambda x: x not in SPLICE_VARS + SV_VARS)
    valid_motif_vars = get_valid_motif_vars(contigs[check_motifs], args)
    contigs.loc[check_motifs, 'valid_motif'] = contigs.variant_id.apply(lambda x: x in valid_motif_vars)

    # check whether exons contain matching novel/partial novel junctions
    spliced_exons = []
    exons = contigs[contigs.variant_type.apply(lambda x: x in NOVEL_BLOCKS)]
    novel_juncs = contigs[contigs.variant_type.apply(lambda x: x in ['PNJ', 'NEJ'])]
    for idx,row in exons.iterrows():
        back_junc = novel_juncs.pos2 == row['pos1']
        front_junc = novel_juncs.pos1 == row['pos2']
        matching_juncs = novel_juncs[np.logical_or(back_junc, front_junc)]
        if len(matching_juncs) > 0:
            spliced_exons.append(row['variant_id'])
    contigs['spliced_exon'] = contigs.variant_id.apply(lambda x: x in spliced_exons)

    # novel exon contigs (spliced, valid motif and large variant sizew)
    is_novel_exon = np.logical_and(contigs.spliced_exon, contigs.valid_motif)
    is_novel_exon = np.logical_and(is_novel_exon, contigs.large_varsize)
    ne_vars = contigs[is_novel_exon].variant_id.values

    # check whether TSVs are within exons
    contigs['overlaps_exon'] = False
    exonic_var = contigs.variant_type.apply(lambda x: x in ['PNJ', 'EE', 'AS'])
    contigs.loc[exonic_var, 'overlaps_exon'] = True
    is_sv = contigs.variant_type.apply(lambda x: x in SV_VARS)
    gene_tree, ex_trees, ex_ref = ac.get_gene_lookup(args.tx_ref_file)
    contigs.loc[is_sv, 'overlaps_exon'] = contigs[is_sv].apply(overlaps_exon, axis=1, args=(ex_trees,))
    sv_vars = contigs[np.logical_and(is_sv, contigs.overlaps_exon)].variant_id.values

    # keep all splice vars
    as_vars = contigs.variant_id.values[contigs.variant_type.apply(lambda x: x in SPLICE_VARS)]

    # ensure retained introns are spliced in some way (to distinguish from pre-mRNAs)
    retained_intron = contigs.variant_type.apply(lambda x: x in ['RI'])
    spliced_ri =  np.logical_and(retained_intron, contigs.is_contig_spliced.values)
    ri_vars = contigs[np.logical_and(spliced_ri, contigs.large_varsize)].variant_id.values

    # collate contigs to keep
    keep_vars = np.unique(np.concatenate([ri_vars, as_vars, ne_vars, sv_vars]))
    contigs['variant_of_interest'] = contigs.variant_id.apply(lambda x: x in keep_vars)
    keep_contigs = contigs[contigs.variant_of_interest].contig_id.values
    contigs = contigs[contigs.contig_id.apply(lambda x: x in keep_contigs)]
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
        exit_with_error(str(exception), EXIT_FILE_IO_ERROR)

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
    keep_contigs = get_contigs_to_keep(args)
    write_output(args, keep_contigs)
    write_bam(args, keep_contigs)

if __name__ == '__main__':
    main()
