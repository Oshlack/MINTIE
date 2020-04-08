import pandas as pd
import re
import sys
import fileinput
from Bio import SeqIO
from argparse import ArgumentParser

def parse_args():
    parser = ArgumentParser(description='Create transcript-to-gene name reference file.')
    parser.add_argument(dest='ref',
                        metavar='REF',
                        type=str,
                        help='''Reference fasta or GTF file. GTF must contain "gene_name" and "transcript_id" fields.
                                Fasta file must contain "gene_symbol:<symbol>" in the header description.''')
    return parser.parse_args()

def make_from_fasta(fasta_file):
    tx2gene = []
    with fileinput.input(fasta_file) as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            tx_id = record.id
            gname = re.search('gene_symbol:([A-Za-z0-9\.\_\-]+)', record.description).group(1)
            tx2gene.append((tx_id, gname))
    tx2gene = pd.DataFrame(tx2gene)
    tx2gene = tx2gene[tx2gene.gene != ''].drop_duplicates()
    tx2gene.to_csv(sys.stdout, index=False, header=False, sep='\t')

def make_from_gtf(gtf):
    tx_ref = pd.read_csv(gtf, comment='#', sep='\t', header=None, low_memory=False)
    if 'transcript' in tx_ref[2].values:
        # makes things faster if the gtf has transcript records
        tx_ref = tx_ref[tx_ref[2] == 'transcript']

    def get_attribute(attributes, attribute_id):
        re_attr = re.search(r'%s "([\w\-\.\/]+)"' % attribute_id, attributes)
        attr = re_attr.group(1) if re_attr else ''
        return attr
    tx_ids = tx_ref[8].apply(lambda x: get_attribute(x, 'transcript_id'))
    genes = tx_ref[8].apply(lambda x: get_attribute(x, 'gene_name'))

    tx2gene = pd.DataFrame.from_dict({'tx': tx_ids, 'gene': genes})
    tx2gene = tx2gene[tx2gene.gene != ''].drop_duplicates()
    tx2gene = tx2gene[['tx', 'gene']]
    tx2gene.to_csv(sys.stdout, index=False, header=False, sep='\t')

def main():
    args = parse_args()
    if args.ref.endswith('.fa') or args.ref.endswith('.fasta'):
        tx2gene = make_from_fasta(args.ref)
    elif args.ref.endswith('.gtf'):
        tx2gene = make_from_gtf(args.ref)
    else:
        print('ERROR: supplied file is not a GTF or fasta file.')
        sys.exit()

if __name__ == '__main__':
    main()
