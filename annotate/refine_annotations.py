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
from argparse import ArgumentParser
from utils import init_logging, exit_with_error

EXIT_FILE_IO_ERROR = 1

MIN_NOVEL_EXON_SIZE = 20
SPLICE_VARS = ['AS', 'PNJ', 'NEJ']
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
    parser.add_argument(dest='contig_out_file',
                        metavar='CONTIG_OUT_FILE',
                        type=str,
                        help='''Contig tsv output file''')
    parser.add_argument(dest='bam_out_file',
                        metavar='BAM_OUT_FILE',
                        type=str,
                        help='''Contig BAM output file''')

    return parser.parse_args()

def get_contigs_to_keep(args):
    '''
    Return contigs matching criteria:
    - novel block size > MIN_NOVEL_EXON_SIZE
    - novel blocks are spliced in some way
    - novel exons have corresponding novel splice sites
    - fusions, SVs and splice variants are not further filtered
    '''
    try:
        cinfo_file = args.contig_info_file
        contigs = pd.read_csv(cinfo_file, sep='\t')
    except IOError as exception:
        exit_with_error(str(exception), EXIT_FILE_IO_ERROR)

    # filter novel blocks: require min size and spliced contig
    is_exon = contigs.variant_type.apply(lambda x: x in NOVEL_BLOCKS)
    is_spliced = contigs.contig_cigar.apply(lambda x: bool(re.search('N|H', x)))
    spliced_exons = np.logical_and(is_spliced, is_exon)
    large_varsize = contigs.contig_varsize > MIN_NOVEL_EXON_SIZE
    novel_exons = contigs[np.logical_and(spliced_exons, large_varsize)]

    # ensure exons contain a matching splice junction
    keep_contigs = []
    novel_juncs = contigs[contigs.variant_type.apply(lambda x: x in ['PNJ', 'NEJ'])]
    for idx,row in novel_exons.iterrows():
        back_junc = novel_juncs.pos2 == row['pos1']
        front_junc = novel_juncs.pos1 == row['pos2']
        matching_juncs = novel_juncs[np.logical_or(back_junc, front_junc)]
        keep_contigs.append(row['contig_id'])

    # keep all SV and splice vars
    svs = contigs[contigs.variant_type.apply(lambda x: x in SV_VARS)]
    splicevars = contigs[contigs.variant_type.apply(lambda x: x in SPLICE_VARS)]

    # ensure retained introns are spliced in some way (to distinguish from pre-mRNAs)
    retained_intron = contigs.variant_type.apply(lambda x: x in ['RI'])
    spliced_ri =  np.logical_and(retained_intron, is_spliced)
    ris = contigs[np.logical_and(spliced_ri, large_varsize)]

    keep_contigs.extend(svs.contig_id)
    keep_contigs.extend(splicevars.contig_id)
    keep_contigs.extend(ris.contig_id)
    keep_contigs = np.unique(keep_contigs)

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
