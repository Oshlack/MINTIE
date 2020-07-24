'''
Module      : get_novel_contigs
Description : Perform contig annotation.
Copyright   : (c) Marek Cmero, Jul 2020
License     : MIT
Maintainer  : github.com/mcmero
Portability : POSIX
Take EC count matrix file, transcriptome
fasta reference and de novo assembled contigs,
and make dummy DE file containing novel contigs.
'''

import pandas as pd
import numpy as np
import os
import sys
from Bio import SeqIO
from argparse import ArgumentParser

def parse_args(args):
    '''
    Parse command line arguments.
    Returns Options object with command line argument values as attributes.
    Will exit the program on a command line error.
    '''
    description = 'Get novel transcripts'
    parser = ArgumentParser(description = description)
    parser.add_argument(dest='ecm_file',
                        metavar='ECM_FILE',
                        type=str,
                        help='''ec_count_matrix.txt file.''')
    parser.add_argument(dest='ref_tx_fasta',
                        metavar='REF_FASTA',
                        type=str,
                        help='''Transcriptome reference fasta.''')
    parser.add_argument(dest='denovo_fasta',
                        metavar='DENOVO_FASTA',
                        type=str,
                        help='''Sample de novo filtered fasta assembly file.''')
    return parser.parse_args(args)

def get_ref_txs(ref_tx_fasta):
    '''
    Get all reference transcript IDs
    from reference transcriptome fasta
    '''
    handle = open(ref_tx_fasta, 'r')
    ref_txs = []
    for record in SeqIO.parse(handle, 'fasta'):
        ref_txs.append(record.id)
    handle.close()
    return ref_txs

def get_tx_ec(ecm, ref_txs, outdir):
    '''
    Make TX > EC lookup table, write the full
    table and return table containing only novel
    contigs
    '''
    # get TX > EC table
    print('Constructing TX > EC table...')
    tx_ec = ecm[['ec_names', 'transcript']].drop_duplicates()
    tx_ec['assembled'] = np.invert(tx_ec.transcript.isin(ref_txs))

    # check whether all TXs in EC are assembled (false if any are reference)
    is_ec_assembled = tx_ec.groupby('ec_names', as_index = False)
    is_ec_assembled = is_ec_assembled.assembled.agg({'all_assembled': all})

    # make full tx_ec table and write this
    tx_ec = tx_ec.merge(is_ec_assembled, on = 'ec_names')
    tx_ec = tx_ec.drop(columns = ['assembled'])
    tx_ec.to_csv('%s/ec_tx_table.txt' % outdir, sep = '\t', index = False)

    # get tx_ec table containing only novel EC txs
    print('Preparing output table...')
    tx_ec = tx_ec[tx_ec.all_assembled]

    # add some extra cols
    contig_count = tx_ec.groupby('ec_names', as_index = False)
    contig_count = contig_count.transcript.agg({'n_contigs_in_ec': len,
                                                'contigs_in_EC': ':'.join})
    tx_ec = tx_ec.merge(contig_count, on = 'ec_names', how = 'left')
    tx_ec = tx_ec.merge(ecm, on = ['ec_names', 'transcript'], how = 'left')
   
    # prepare for output
    sample_name = ecm.columns.values[0]
    tx_ec = tx_ec.rename({'transcript': 'contig', sample_name: 'case_reads'}, axis=1)
    tx_ec = tx_ec[['ec_names', 'contig', 'n_contigs_in_ec', 'contigs_in_EC', 'case_reads']]

    return tx_ec

def main():
    args = parse_args(sys.argv[1:])
    try:
        print('Reading in EC count matrix file...')
        ecm = pd.read_csv(args.ecm_file, sep = '\t')
        print('Fetching reference transcripts...')
        ref_txs = get_ref_txs(args.ref_tx_fasta)
    except IOError as message:
        print("{} ERROR: {}, exiting".format("get_novel_contigs", message), file=sys.stderr)
        sys.exit(1)

    # get tx_ec table
    outdir = os.path.dirname(args.ecm_file)
    tx_ec = get_tx_ec(ecm, ref_txs, outdir)
    
    # write dummy DE file
    tx_ec.to_csv('%s/eq_classes_de.txt' % outdir, sep = '\t', index = False)

if __name__ == '__main__':
    main()
