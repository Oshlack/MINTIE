#########################################################
# Author: Marek Cmero
#########################################################

import argparse
import pandas as pd
import numpy as np
import gffutils
from Bio import SeqIO

import ipdb

parser = argparse.ArgumentParser()
parser.add_argument(dest='contigs_fasta',
                    help='''Fasta file containing contig sequences''')
parser.add_argument(dest='novel_contigs_file',
                    help='''Annotated novel contigs file''')
parser.add_argument(dest='block_bed',
                    help='''Sequence blocks used to construct transcript sequences.''')
parser.add_argument(dest='block_fasta',
                    help='''Fasta file containing sequences (contructed from the blocks bed file)''')
parser.add_argument(dest='output_fasta',
                    help='''Output supertranscript file''')

args        = parser.parse_args()
con_fasta   = args.contigs_fasta
nc_file     = args.novel_contigs_file
blocks      = args.block_bed
block_fasta = args.block_fasta
st_file     = args.output_fasta

novel_contigs = pd.read_csv(nc_file, sep='\t')
#novel_contigs = novel_contigs[novel_contigs.variant != 'soft-clip']
exon_df =  pd.read_csv(blocks, sep='\t', header = None, names = ['chrom', 'start', 'end', 'name', 'value', 'strand'])
exon_df[['gene', 'blocks']] = exon_df['name'].str.split(' ', expand=True)

block_seqs = {}
handle = open(block_fasta, 'r')
for record in SeqIO.parse(handle, 'fasta'):
    block_seqs[record.id] = str(record.seq)
handle.close()

contig_seqs = {}
handle = open(con_fasta, 'r')
for record in SeqIO.parse(handle, 'fasta'):
    contig_seqs[record.id] = record.seq
handle.close()

gene_out = pd.DataFrame()
genes = np.unique([gn for gene in novel_contigs.gene.values for gn in gene.split('|')])
for gene in genes:
    nc = novel_contigs[novel_contigs.gene == gene]
    gene_df = exon_df[exon_df.gene == gene]
    gene_out = gene_df.copy()
    if len(gene_out) == 0:
        continue
    antisense = True if gene_out.strand.values[0]=='-' else False

    contigs_with_novel_bits = nc[nc.contig_varsize.values > 0]
    if any(nc.contig_varsize.values > 0):
        for idx, contig_row in contigs_with_novel_bits.iterrows():
            seq = contig_seqs[contig_row['contig']]
            csize = contig_row.contig_varsize
            start = (len(seq) - contig_row.contig_pos2) if antisense else contig_row.contig_pos1
            end = start + csize

            gpos1, gpos2 = contig_row.genome_pos1, contig_row.genome_pos2
            blocks_affected = gene_df[np.logical_and(gene_df.start < gpos1, gene_df.end > gpos2)]

            chrom = 'MT' if contig_row.chrom1 == 'chrM' else contig_row.chrom1.split('chr')[1]
            novel_seq = seq[start:end]
            novel_seq_info = (chrom, gpos1, gpos2, gene_out.strand.values[0])
            block_seqs['%s:%d-%d(%s)' % novel_seq_info] = novel_seq

            assert len(blocks_affected) <= 1
            for idx, block in blocks_affected.iterrows():
                seq = block_seqs['%s:%d-%d(%s)' % (chrom, block.start, block.end, block.strand)]
                left_seq = seq[:block.end-gpos2] if antisense else seq[:gpos1-block.start]
                right_seq = seq[block.end-gpos1:] if antisense else seq[gpos2-block.start:]

                novel_block = pd.DataFrame([[block.chrom, gpos1, gpos2, '%s novel' % gene,
                                            block.value, block.strand, block.gene, 'novel']], columns=gene_out.columns)
                add_blocks = [novel_block]

                if gpos1 - block.start != 0:
                    block_seqs['%s:%d-%d(%s)' % (chrom, block.start, gpos1, block.strand)] = right_seq if antisense else left_seq
                    left_block = pd.DataFrame([[block.chrom, block.start, gpos1, block['name'],
                                                block.value, block.strand, block.gene, block['blocks']]], columns=gene_out.columns)
                    add_blocks.append(left_block)

                if block.end - gpos2 != 0:
                    block_seqs['%s:%d-%d(%s)' % (chrom, gpos2, block.end, block.strand)] = left_seq if antisense else right_seq
                    right_block = pd.DataFrame([[block.chrom, gpos2, block.end, block['name'],
                                                 block.value, block.strand, block.gene, block['blocks']]], columns=gene_out.columns)
                    add_blocks.append(right_block)

                gene_out = gene_out.append(add_blocks, ignore_index=True)
                gene_out = gene_out[gene_out.index!=block.name]

    gene_out = gene_out.sort_values(by=['start', 'end'], ascending=False) if antisense else gene_out.sort_values(by=['start','end'])
    seqs, names = [], []
    for idx,x in gene_out.iterrows():
        names.append(x['blocks'])
        seq = str(block_seqs['%s:%d-%d(%s)' % (x.chrom, x.start, x.end, x.strand)])
        seqs.append(seq)

    seq_ends = np.cumsum([len(s) for s in seqs])
    seq_starts = np.concatenate([[0], seq_ends[:-1]])
    segs = ['%s-%s' % (s1+1, s2) for s1,s2 in zip(seq_starts, seq_ends)]

    header = '>%s segs:%s names:%s\n' % (gene, ','.join(segs), ','.join(names))
    sequence = ''.join(seqs) + '\n'
    with open(st_file, 'a') as st_fasta:
        st_fasta.writelines([header, sequence])
