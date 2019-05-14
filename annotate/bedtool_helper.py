'''
Module      : bedtool_helper
Description : Helper functions for pybedtools
Copyright   : (c) Marek Cmero, Mar 2019
License     : TBD
Maintainer  : MAREK.CMERO@MCRI.EDU.AU
Portability : POSIX
'''
import tempfile
import pandas as pd
from Bio import SeqIO
from pybedtools import BedTool
from pybedtools.featurefuncs import extend_fields

def featuretype_filter(feature, featuretype):
    '''
    from http://daler.github.io/pybedtools/3-brief-examples.html
    Only passes features with the specified *featuretype*
    '''
    if feature[2] == featuretype:
        return True
    return False

def subset_featuretypes(g, featuretype):
    '''
    from http://daler.github.io/pybedtools/3-brief-examples.html
    Returns the filename containing only `featuretype` features.
    '''
    return g.filter(featuretype_filter, featuretype).saveas().fn

def add_strand(exon, strand):
    '''
    add strand to exon block
    '''
    exon = extend_fields(exon, 6)
    exon.strand = strand
    return exon

def get_block_seqs(exons):
    '''
    get sequences from exon blocks and
    return block sequences dictionary
    '''
    block_seqs = {}
    with tempfile.NamedTemporaryFile() as fa_tmp:
        fa_tmp.write(bytes(open(exons.seqfn).read(), 'utf-8'))
        fa_tmp.flush()

        for record in SeqIO.parse(fa_tmp.name, 'fasta'):
            block_seqs[record.id] = str(record.seq)

    return(block_seqs)

def get_merged_exons(genes, gtf, genome_fasta, strand):
    '''
    get all exons from specified genes, merging any
    overlapping exonic regions, also return their
    respective sequences in a dictionary object
    '''
    gene_gtf = gtf[gtf.gene.isin(genes)]
    if len(gene_gtf) == 0:
        return pd.DataFrame(), {}
    gene_gtf = gene_gtf.drop('gene', axis=1)
    gene_strand = gene_gtf.strand.values[0]
    strand = gene_strand if strand == '' else strand

    with tempfile.NamedTemporaryFile(mode='r+') as temp_gtf:
        gene_gtf.to_csv(temp_gtf.name, index=False, header=False, sep='\t')

        # load gene GTF info, extract and merge exons
        g = BedTool(temp_gtf.name)
        exons = BedTool(subset_featuretypes(g, 'exon'))
        exons = exons.remove_invalid().sort().merge()

        exseq = exons.each(add_strand, strand)
        exseq = exseq.sequence(fi=genome_fasta, s=True)
        block_seqs = get_block_seqs(exseq)

        blocks = pd.DataFrame()
        with tempfile.NamedTemporaryFile(mode='r+') as temp_exons:
            exons.saveas(temp_exons.name)
            blocks = pd.read_csv(temp_exons, header=None, sep='\t', names=['chr', 'start', 'end'])

        if type(genes) == str:
            blocks['name'] = genes
        else:
            blocks['name'] = '|'.join(list(genes))

        blocks['score'] = '.'
        blocks['strand'] = strand

    blocks['chr'] = blocks['chr'].map(str)
    blocks.start = blocks.start.map(int)
    blocks.end = blocks.end.map(int)

    block_names = []
    # reverse numbering if the gene is on the reverse strand
    if gene_strand == '-':
        block_names = ['|' + str(i) for i in reversed(range(1, len(blocks)+1))]
    else:
        block_names = ['|' + str(i) for i in range(1, len(blocks)+1)]
    blocks['name'] = blocks['name'] + block_names

    return(blocks, block_seqs)
