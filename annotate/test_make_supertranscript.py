import pytest
import pandas as pd
import numpy as np
import make_supertranscript as ms


@pytest.mark.parametrize('seq,expected', [('AGTC', 'GACT'),
                                          ('NAGTC', 'GACTN'),
                                          ('A', 'T'),
                                          ('C', 'G'),
                                          ('T', 'A'),
                                          ('G', 'C')])
def test_reverse_complement(seq, expected):
    assert ms.reverse_complement(seq) == expected

def test_get_contig_genes():
    con_info = {'overlapping_genes': ['A:B', 'B:', 'C']}
    con_info = pd.DataFrame.from_dict(con_info)
    assert list(ms.get_contig_genes(con_info)) == ['A', 'B']

def test_get_contig_strand():
    contigs = {'pos1': ['chr1:100(+)', 'chr1:200(-)', 'chr1:400(+)'],
               'pos2': ['chr1:100(-)', 'chr1:200(+)', 'chr1:500(+)'],
               'variant_id': ['A', 'B', 'C'],
               'partner_id': ['B', 'A', '']}
    contigs = pd.DataFrame.from_dict(contigs)
    assert ms.get_contig_strand(contigs, 'A') == '+'
    assert ms.get_contig_strand(contigs, 'B') == '-'
    assert ms.get_contig_strand(contigs, 'C') == '+'
    assert ms.get_contig_strand(contigs, 'X') == '.'
