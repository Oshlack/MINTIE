import pytest
import pandas as pd
import annotate_contigs as ac
from intervaltree import Interval, IntervalTree

# dummy global args
args = type('argparse.Namespace', (object,),
            {'minClip': 20, 'minGap': 7, 'minMatch': '30, 0.3'})()
ac.set_globals(args)

# dummy read class
class Read:
    reads = []
    def __init__(self, reference_name, blocks, query_name):
        self.reference_name = reference_name
        self.blocks = blocks
        self.query_name = query_name
        Read.reads.append(self)
    def get_blocks(self):
        return self.blocks
    def find(self, query_name):
        query_reads = [read for read in Read.reads if read.query_name == query_name]
        return query_reads

@pytest.mark.parametrize('letter,next_letter', [('a', 'b'),
                                                ('g', 'h'),
                                                ('z', 'A'),
                                                ('Z', 'Za'),
                                                ('ZZ', 'ZZa'),
                                                ('Za', 'Zb')])
def test_get_next_letter(letter, next_letter):
    assert ac.get_next_letter(letter) == next_letter

def test_get_next_id():
    assert ac.get_next_id('k49_123') == 'k49_123a'
    assert ac.get_next_id('k49_123') == 'k49_123b'

def test_get_juncs():
    tx_in = {'chrom': '1',
             'exonStarts': '100,300,400',
             'exonEnds': '200,350,600'}
    result = [('1', '200', '300'),
              ('1', '350', '400')]
    assert ac.get_juncs(tx_in) == result

def test_get_tx_juncs():
    read = Read('chr1',
                [(100, 200), (200, 300),
                 (305, 400), (500, 600)],
                 'A')
    assert ac.get_tx_juncs(read) == [('chr1', 400, 500)]

@pytest.mark.parametrize('chrom,expected', [('chr1', 'intervals_a'),
                                            ('chr2', 'intervals_b'),
                                            ('chrM', 'intervals_c'),
                                            ('MT', 'intervals_c'),
                                            ('chr3', 'intervals_d'),
                                            ('3', 'intervals_d'),
                                            ('1', 'intervals_a'),
                                            ('GL123', None)])
def test_get_chrom_ref_tree(chrom, expected):
    ref_trees = {'chr1': 'intervals_a',
                 'chr2': 'intervals_b',
                 'chrM': 'intervals_c',
                 '3': 'intervals_d'}
    assert ac.get_chrom_ref_tree(chrom, ref_trees) == expected


# note that for this test, reads are cumulatively added to the Read object
@pytest.mark.parametrize('read,expected', [('chr1:90-105:A', True),
                                           ('chr1:250-290:A', True),
                                           ('chr1:250-290:B', False),
                                           ('chrX:100-200:C', False)])
def test_do_any_read_blocks_overlap_exons(read, expected):
    ex_trees, ref_tree = {}, IntervalTree()
    coords = [(100, 200),
              (300, 400)]
    for s,e in coords:
        ref_tree.addi(s, e)
    ex_trees['chr1'] = ref_tree

    chrom, coords, name = read.split(':')
    start, end = int(coords.split('-')[0]), int(coords.split('-')[1])
    read = Read(chrom, [(start, end)], name)

    assert ac.do_any_read_blocks_overlap_exons(read, ex_trees, read) == expected


@pytest.mark.parametrize('read,expected', [('chr1:90-105:A', 'GeneA'),
                                           ('chr1:250-290:A', ''),
                                           ('chr1:350-389:B', 'GeneB'),
                                           ('chr1:350-450:B', 'GeneB|GeneC'),
                                           ('chrX:100-200:C', '')])
def test_get_overlapping_genes(read, expected):
    gn_trees, ref_tree = {}, IntervalTree()
    coords = [(100, 200, 'GeneA'),
              (300, 400, 'GeneB'),
              (390, 500, 'GeneC')]
    for s,e,g in coords:
        ref_tree.addi(s, e, g)
    gn_trees['chr1'] = ref_tree

    chrom, coords, name = read.split(':')
    start, end = int(coords.split('-')[0]), int(coords.split('-')[1])
    read = Read(chrom, [(start, end)], name)

    assert ac.get_overlapping_genes(read, gn_trees) == expected
