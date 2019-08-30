import pytest
import pandas as pd
import annotate_contigs as ac

# dummy global args
args = type('argparse.Namespace', (object,),
            {'minClip': 20, 'minGap': 7, 'minMatch': '30, 0.3'})()
ac.set_globals(args)

# dummy read class
class Read:
    def __init__(self, reference_name, blocks):
        self.reference_name = reference_name
        self.blocks = blocks
    def get_blocks(self):
        return self.blocks

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
                 (305, 400), (500, 600)])
    assert ac.get_tx_juncs(read) == [('chr1', 400, 500)]

@pytest.mark.parametrize('chrom,expected', [('chr1', 'intervals_a'),
                                            ('chr2', 'intervals_b'),
                                            ('chrM', 'intervals_c'),
                                            ('MT', 'intervals_c'),
                                            ('chr3', 'intervals_d'),
                                            ('3', 'intervals_d'),
                                            ('1', 'intervals_a')])
def test_get_chrom_ref_tree(chrom, expected):
    ref_trees = {'chr1': 'intervals_a',
                 'chr2': 'intervals_b',
                 'chrM': 'intervals_c',
                 '3': 'intervals_d'}
    assert ac.get_chrom_ref_tree(chrom, ref_trees) == expected
