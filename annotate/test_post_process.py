import pytest
import pandas as pd
import numpy as np
import post_process as pp

# dummy classes
class Read:
    def __init__(self, contig, start, end, name):
        self.contig = contig
        self.reference_start = start
        self.reference_end = end
        self.query_name = name
class AlignmentFile:
    def __init__(self, reads):
        self.reads = []
        for read in reads:
            contig, start, end, name = read
            self.reads.append(Read(contig, start, end, name))
    def fetch(self, contig, start, end, until_eof=True):
        # dummy method
        return self.reads

test_reads = [('1', 100, 200, 'A'),
              ('1', 105, 205, 'B'),
              ('1', 200, 300, 'C')]

def test_parse_args():
    args = pp.parse_args(['sample',
                         'contig_info.tsv',
                         'de_results.tsv',
                         'estimated_VAF.txt',
                         '--gene_filter', 'gene_filter.txt',
                         '--var_filter', 'FUS INS DEL'])
    assert args.sample == 'sample'
    assert args.contig_info == 'contig_info.tsv'
    assert args.gene_filter == 'gene_filter.txt'
    assert args.var_filter == ['FUS INS DEL']

@pytest.mark.parametrize('gene,expected', [('A', ['A']),
                                            ('A|B', ['A', 'B']),
                                            ('A:B', ['A', 'B']),
                                            ('A:B|C', ['A', 'B', 'C'])])
def test_get_all_genes(gene, expected):
    pp.get_all_genes(gene) == expected

def test_filter_by_gene():
    gene_filter = {0: ['B', 'C', 'E', 'X']}
    gene_filter = pd.DataFrame.from_dict(gene_filter)

    contigs = {'overlapping_genes': ['A', 'A:B', 'C|E', 'Y', 'X']}
    contigs = pd.DataFrame.from_dict(contigs)

    result = pp.filter_by_gene(contigs, gene_filter).overlapping_genes.values
    assert list(result) == ['A:B', 'C|E', 'X']

def test_add_de_info():
    de_results = {'contig': ['k49_123', 'k79_234', 'k49_789'],
                  'contigs': ['k49_123:k79_234', 'k49_123:k79_234', 'k49_789'],
                  'case_reads': [200, 200, 50],
                  'controls_total_reads': [5, 5, 10],
                  'logFC': [5, 5, 2],
                  'F': [50, 50, 25],
                  'logCPM': [10, 10, 5],
                  'n_contigs_in_ec': [2, 2, 1],
                  'PValue': [0.001, 0.001, 0.003],
                  'FDR': [0.01, 0.01, 0.03],
                  'ec_names': ['ec1', 'ec1', 'ec3']}
    de_results = pd.DataFrame.from_dict(de_results)
    contigs = {'contig_id': ['k49_123', 'k79_234', 'k49_789']}
    contigs = pd.DataFrame.from_dict(contigs)

    results = {'contig_id': ['k49_123', 'k79_234', 'k49_789'],
               'contigs_in_EC': ['k49_123:k79_234', 'k49_123:k79_234', 'k49_789'],
               'case_reads': [200, 200, 50],
               'controls_total_reads': [5, 5, 10],
               'logFC': [5, 5, 2],
               'F': [50, 50, 25],
               'logCPM': [10, 10, 5],
               'n_contigs_in_ec': [2, 2, 1],
               'PValue': [0.001, 0.001, 0.003],
               'FDR': [0.01, 0.01, 0.03],
               'ec_names': ['ec1', 'ec1', 'ec3']}
    results = pd.DataFrame.from_dict(results)

    assert all(pp.add_de_info(contigs, de_results) == results)

@pytest.mark.parametrize('gene,expected', [('A', 'A'),
                                           ('A|B', 'A'),
                                           ('A:B', 'A|B'),
                                           ('A|B|C:X|Y|Z', 'A|X')])
def test_get_short_gene_name(gene, expected):
    assert pp.get_short_gene_name(gene) == expected

def test_reformat_fields():
    contigs = { 'contig_id': ['1', '2'],
                'variant_id': ['1a', '2a'],
                'partner_id': ['', ''],
                'vars_in_contig': [1, 1],
                'pos1': ['chr1:100(+)', 'chr2:200(-)'],
                'pos2': ['chr1:200(+)', 'chr2:250(-)'],
                'varsize': [100, 0],
                'cpos': [0, 0],
                'contig_varsize': [0, 0],
                'contig_len': [200, 150],
                'contig_cigar': ['100M50N100M', '50M50N50M'],
                'variant_type': ['NE', 'NEJ'],
                'overlapping_genes': ['A', 'B'],
                'large_varsize': [True, False],
                'is_contig_spliced': [True, True],
                'spliced_exon': [True, False],
                'overlaps_exon': [False, False],
                'overlaps_gene': [True, True],
                'valid_motif': [False, True],
                'sample': ['S1', 'S1'],
                'logFC': [10, 10],
                'logCPM': [10, 9],
                'case_CPM': [15, 12],
                'PValue': [0.001, 0.002],
                'FDR': [0.01, 0.02],
                'ec_names': ['ec1', 'ec2'],
                'n_contigs_in_ec': [1, 1],
                'contigs_in_EC': ['1', '2'],
                'case_reads': [100, 200],
                'controls_total_reads': [0, 5],
                'TPM': [1000, 2000],
                'mean_WT_TPM': [1400, 1400],
                'VAF': [1000./2400, 2000./3400],
                'unique_contig_ID': ['S1|1|A', 'S1|2|B'] }
    contigs = pd.DataFrame.from_dict(contigs)
    contigs = pp.reformat_fields(contigs)

    assert list(contigs['chr1'].values) == ['chr2', 'chr1']
    assert list(contigs['chr2'].values) == ['chr2', 'chr1']
    assert list(contigs['pos1'].values) == [200, 100]
    assert list(contigs['pos2'].values) == [250, 200]
    assert list(contigs['strand1'].values) == ['-', '+']
    assert list(contigs['strand2'].values) == ['-', '+']
