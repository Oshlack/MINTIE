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
    contigs = {'pos1': ['chr1:100(+)', 'chr1:400(+)'],
               'pos2': ['chr1:200(-)', 'chr1:500(+)'],
               'variant_id': ['A', 'C'],
               'partner_id': ['B', '']}
    contigs = pd.DataFrame.from_dict(contigs)
    assert ms.get_contig_strand(contigs, 'A') == '+'
    assert ms.get_contig_strand(contigs, 'B') == '-'
    assert ms.get_contig_strand(contigs, 'C') == '+'
    assert ms.get_contig_strand(contigs, 'X') == '.'

def test_get_gene_strands():
    gtf = {'gene': ['A', 'B', 'C', 'D', ''],
           'strand': ['+', '+', '-', '', '+']}
    gtf = pd.DataFrame.from_dict(gtf)
    assert ms.get_gene_strands(gtf, ['A', 'B', 'C', 'D']) == ['+', '+', '-', '']

@pytest.mark.parametrize('strands,expected', [((['+', '-'], ['+', '-']), ['+', '-']),
                                              ((['-', '-'], ['+', '+']), ['+', '+']),
                                              ((['+', '-'], ['-', '-']), ['+', '-']),
                                              ((['-', '-'], ['-', '+']), ['-', '-'])])
def test_get_strand_info(strands, expected):
    cstrands, gstrands = strands
    pos1, pos2 = 'chr1:100(%s)' % cstrands[0], 'chr1:200(%s)' % cstrands[1]
    contigs = {'pos1': [pos1],
               'pos2': [pos2],
               'variant_id': ['A'],
               'partner_id': ['B'],
               'variant_type': ['FUS']}
    contigs = pd.DataFrame.from_dict(contigs)
    assert ms.get_strand_info(contigs, gstrands) == expected

def test_get_output_files():
    sample, outdir = 'A', '.'
    genome_bed, st_block_bed, st_gene_bed, st_fasta = ms.get_output_files(sample, outdir)
    assert genome_bed == './A_genome.bed'
    assert st_block_bed == './A_blocks_supertranscript.bed'
    assert st_gene_bed == './A_genes_supertranscript.bed'
    assert st_fasta == './A_supertranscript.fasta'
