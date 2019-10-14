import pytest
import pandas as pd
import numpy as np
import post_process as pp
import count_junction_reads as cjr

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
                         'st_bed.bed',
                         'cont_align.bam',
                         'read_align.bam',
                         '--gene_filter', 'gene_filter.txt',
                         '--var_filter', 'FUS INS DEL'])
    assert args.sample == 'sample'
    assert args.contig_info == 'contig_info.tsv'
    assert args.st_bed == 'st_bed.bed'
    assert args.cont_align == 'cont_align.bam'
    assert args.read_align == 'read_align.bam'
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
    de_results = {'contig': ['A', 'A', 'B'],
                  'contigs': ['k49_123', 'k79_234', 'k49_789'],
                  'case_reads': [100, 200, 50],
                  'controls_total_reads': [5, 20, 10],
                  'logFC': [5, 6, 2],
                  'F': [45, 50, 25],
                  'logCPM': [10, 15, 5],
                  'n_contigs_in_ec': [1, 2, 1],
                  'PValue': [0.005, 0.001, 0.003],
                  'FDR': [0.05, 0.01, 0.03],
                  'ec_names': ['ec1', 'ec2', 'ec3']}
    de_results = pd.DataFrame.from_dict(de_results)
    contigs = {'contig_id': ['A', 'B']}
    contigs = pd.DataFrame.from_dict(contigs)

    results = {'contig_id': ['A', 'B'],
               'contigs_in_EC': ['k49_123,k79_234', 'k49_789'],
               'case_reads': [300, 50],
               'controls_total_reads': [25, 10],
               'logFC': [6, 2],
               'F': [50, 25],
               'logCPM': [15, 5],
               'n_contigs_in_ec': [2, 1],
               'PValue': [0.001, 0.003],
               'FDR': [0.01, 0.03],
               'ec_names': ['ec1,ec2', 'ec3']}
    results = pd.DataFrame.from_dict(results)

    assert all(pp.add_de_info(contigs, de_results) == results)

@pytest.mark.parametrize('gene,expected', [('A', 'A'),
                                           ('A|B', 'A'),
                                           ('A:B', 'A|B'),
                                           ('A|B|C:X|Y|Z', 'A|X')])
def test_get_short_gene_name(gene, expected):
    assert pp.get_short_gene_name(gene) == expected

def test_make_junctions():
    st_blocks = {'start': [100, 200, 300],
                 'end': [150, 202, 311]}
    st_blocks = pd.DataFrame.from_dict(st_blocks)
    result = {'start': [100, 150, 200, 300, 311],
              'end': [100, 150, 202, 300, 311]}
    result = pd.DataFrame.from_dict(result)
    juncs = pp.make_junctions(st_blocks).sort_values('start').reset_index(drop=True)
    assert all(result == juncs)

def test_get_crossing_reads():
    bamf = AlignmentFile(test_reads)
    assert list(cjr.get_crossing_reads('1', 110, 117, bamf)) == ['A', 'B']

@pytest.mark.parametrize('junc,expected', [(('1', 110, 116), 2),
                                           (('1', 250, 250), 1)])
def test_get_read_counts(junc, expected):
    bamf = AlignmentFile(test_reads)
    contig, start, end = junc
    junc = {'contig': [contig],
             'start': [start],
             'end': [end]}
    junc = pd.DataFrame.from_dict(junc)
    assert cjr.get_read_counts(bamf, junc).crossing.values[0] == expected

def test_get_read_support():
    st_bed = {'contig': ['s1|A', 's1|A'],
              'start': [110, 250],
              'end': [116, 250]}
    st_bed = pd.DataFrame.from_dict(st_bed)
    contigs = {'contig_id': ['A'],
               'sample': ['s1'],
               'contig': ['1']}
    contigs = pd.DataFrame.from_dict(contigs)
    bamf = AlignmentFile(test_reads)
    assert pp.get_read_support(contigs, bamf, st_bed).crossing_reads.values[0] == '2,1'
