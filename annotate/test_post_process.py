import pytest
import post_process as pp

@pytest.mark.parametrize('input,expected', [('A', ['A']),
                                            ('A|B', ['A', 'B']),
                                            ('A:B', ['A', 'B']),
                                            ('A:B|C', ['A', 'B', 'C'])])
def test_get_all_genes(input, expected):
    pp.get_all_genes(input) == expected
