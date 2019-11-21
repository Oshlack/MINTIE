import pandas as pd
import re
import sys
from argparse import ArgumentParser

def parse_args():
    parser = ArgumentParser(description='Create transcript-to-gene name reference file.')
    parser.add_argument(dest='gtf',
                        metavar='GTF',
                        type=str,
                        help='''GTF file. Must contain "gene_name" and "transcript_id" fields.''')
    return parser.parse_args()

def get_attribute(attributes, attribute_id):
    re_attr = re.search(r'%s "([\w\-\.\/]+)"' % attribute_id, attributes)
    attr = re_attr.group(1) if re_attr else ''
    return attr

def main():
    args = parse_args()
    tx_ref = pd.read_csv(args.gtf, comment='#', sep='\t', header=None, low_memory=False)
    if 'transcript' in tx_ref[2].values:
        # makes things faster if the gtf has transcript records
        tx_ref = tx_ref[tx_ref[2] == 'transcript']

    tx_ids = tx_ref[8].apply(lambda x: get_attribute(x, 'transcript_id'))
    genes = tx_ref[8].apply(lambda x: get_attribute(x, 'gene_name'))

    tx2gene = pd.DataFrame.from_dict({'tx': tx_ids, 'gene': genes})
    tx2gene = tx2gene[tx2gene.gene != ''].drop_duplicates()
    tx2gene.to_csv(sys.stdout, index=False, header=False, sep='\t')

if __name__ == '__main__':
    main()
