#########################################################
# Author: Marek Cmero
#########################################################

import argparse
import pandas as pd
import numpy as np
import gffutils
import math
import os
from Bio import SeqIO

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
sample      = os.path.dirname(con_fasta).split('/')[-1].split('_')[0]

# for reverse-complementing
lookup = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

def reverse_complement(seq):
    if seq == '':
        return ''
    if type(seq) == float and math.isnan(seq):
        return ''
    seq = seq[::-1]
    seq = ''.join([lookup[base] for base in list(seq)])
    return(seq)

def split_block(nc_row, gene_out, seq, block, gpos1, gpos2):
    seq = block_seqs['%s:%d-%d(%s)' % (chrom, block.start, block.end, block.strand)]
    left_seq = seq[:block.end-gpos2] if antisense else seq[:gpos1-block.start]
    right_seq = seq[block.end-gpos1:] if antisense else seq[gpos2-block.start:]
    add_blocks = []
    variant = nc_row.variant.replace(' ', '_')
    if nc_row.contig_varsize > 0:
        novel_block = pd.DataFrame([[block.chrom, gpos1, gpos2, '%s %s' % (gene, variant),
                                    block.value, block.strand, block.gene, variant]], columns=gene_out.columns)
        add_blocks = [novel_block]
    else:
        var_seq = seq[block.end-gpos2:block.end-gpos1] if antisense else seq[gpos1-block.start:gpos2-block.start] #actually a reference sequence
        block_seqs['%s:%d-%d(%s)' % novel_seq_info] = str(var_seq)
        var_block = pd.DataFrame([[block.chrom, gpos1, gpos2, '%s %s|%s' % (gene, variant, block['blocks']),
                                   block.value, block.strand, block.gene, '%s|%s' % (variant, block['blocks'])]], columns=gene_out.columns)
        add_blocks = [var_block]

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

    gene_out = gene_out[gene_out.index!=block.name]
    gene_out = gene_out.append(add_blocks, ignore_index=True)

    return gene_out, block_seqs

def write_gene(gene_out, block_seqs, name, st_file):
    seqs, names = [], []
    for idx,x in gene_out.iterrows():
        names.append(x['blocks'])
        seq = str(block_seqs['%s:%d-%d(%s)' % (x.chrom, x.start, x.end, x.strand)])
        seqs.append(seq)

    seg_ends = np.cumsum([len(s) for s in seqs])
    seg_starts = np.concatenate([[0], seg_ends[:-1]])
    segs = ['%s-%s' % (s1+1, s2) for s1,s2 in zip(seg_starts, seg_ends)]

    header = '>%s segs:%s names:%s\n' % (name, ','.join(segs), ','.join(names))
    sequence = ''.join(seqs) + '\n'
    with open(st_file, 'a') as st_fasta:
        st_fasta.writelines([header, sequence])

novel_contigs = pd.read_csv(nc_file, sep='\t')
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

tmp = novel_contigs.gene.str.split('|').apply(pd.Series, 1).stack()
tmp.index = tmp.index.droplevel(-1)
tmp = pd.DataFrame(tmp, columns=['gene_name'])
tmp.name = 'gene_name'
novel_contigs = novel_contigs.join(tmp)
del novel_contigs['gene_name']

gene_out = pd.DataFrame()
genes = np.unique([gn for gene in novel_contigs.gene.values for gn in gene.split('|')])
for gene in genes:
    nc = novel_contigs[[gene in gn.split('|') for gn in novel_contigs.gene.values]]
    gene_df = exon_df[[gene in gn.split('|') for gn in exon_df.gene.values]]
    gene_out = gene_df.copy()
    if len(gene_out) == 0:
        continue

    antisense = True if gene_out.strand.values[0]=='-' else False
    gene_df = gene_df.drop_duplicates().sort_values(by=['start', 'end'], ascending=False).reset_index(drop=True) \
                    if antisense else gene_df.drop_duplicates().sort_values(by=['start','end']).reset_index(drop=True)

    var_count = 0
    for idx, nc_row in nc.iterrows():
        # not necessary to annotate deletions or junctions with no novel sequence
        if nc_row.variant == 'deletion':
            continue
        if nc_row.variant in ['novel junction', 'fusion']:
            if type(nc_row.variant_seq) is str and len(nc_row.variant_seq) == 0:
                continue
            elif type(nc_row.variant_seq) is float and math.isnan(nc_row.variant_seq):
                continue

        seq = contig_seqs[nc_row['contig']]
        csize = nc_row.contig_varsize
        start = (len(seq) - nc_row.contig_pos2) if antisense else nc_row.contig_pos1
        end = start + csize

        gpos1, gpos2 = nc_row.genome_pos1, nc_row.genome_pos2
        blocks_affected = []

        chrom = 'MT' if nc_row.chrom1 == 'chrM' else nc_row.chrom1.split('chr')[1]
        gene_strand = '-' if antisense else '+'
        novel_seq = nc_row.variant_seq

        if nc_row.variant not in ['novel junction', 'fusion'] and antisense:
            novel_seq = reverse_complement(nc_row.variant_seq)
        novel_seq_info = (chrom, gpos1, gpos2, gene_out.strand.values[0])

        blocks_affected = pd.DataFrame()
        if nc_row.variant == 'fusion':
            chrom1 = 'MT' if nc_row.chrom1 == 'chrM' else nc_row.chrom1.split('chr')[1]
            chrom2 = 'MT' if nc_row.chrom2 == 'chrM' else nc_row.chrom2.split('chr')[1]

            novel_seq_info = (chrom1, gpos1, gpos1, gene_out.strand.values[0])
            block_seqs['%s:%d-%d(%s)' % novel_seq_info] = novel_seq

            novel_seq_info = (chrom2, gpos2, gpos2, gene_out.strand.values[0])
            block_seqs['%s:%d-%d(%s)' % novel_seq_info] = novel_seq

            if chrom1 == chrom2:
                ba1 = gene_df[np.logical_and(gene_df.start < gpos1, gene_df.end > gpos1)]
                if len(ba1) > 0:
                    block = gene_df.loc[ba1.index.values[0]]
                    gene_out, block_seqs = split_block(nc_row, gene_out, seq, block, gpos1, gpos1)

                ba2 = gene_out[np.logical_and(gene_out.start < gpos2, gene_out.end > gpos2)]
                if len(ba2) > 0:
                    block = gene_out.loc[ba2.index.values[0]]
                    gene_out, block_seqs = split_block(nc_row, gene_out, seq, block, gpos2, gpos2)
            else:
                chrom = chrom1 if chrom1 == gene_df.chrom.values[0] else chrom2
                gpos = gpos1 if chrom1 == gene_df.chrom.values[0] else gpos2

                blocks_affected = gene_df[np.logical_and(gene_df.start < gpos, gene_df.end > gpos)]
                if len(blocks_affected) > 0:
                    block = blocks_affected.loc[blocks_affected.index.values[0]]
                    gene_out, block_seqs = split_block(nc_row, gene_out, seq, block, gpos, gpos)
        elif nc_row.variant == 'novel junction':
            # determine whether novel sequence should be inserted on left or right of an exon
            on_sense = nc_row.contig_align_strand == '+'
            seq_on_contig_start = (nc_row.contig_pos1 == 0)
            left_side = (on_sense and seq_on_contig_start) or (not on_sense and not seq_on_contig_start)

            left_affected_block = gene_df[np.logical_and(gene_df.start < gpos1, gene_df.end > gpos1)]
            right_affected_block = gene_df[np.logical_and(gene_df.start < gpos2, gene_df.end > gpos2)]

            if len(left_affected_block) > 0 or len(right_affected_block) > 0:
                affected_blocks = left_affected_block.append(right_affected_block).drop_duplicates()
                if left_side:
                    gpos2 = affected_blocks.start.values[0]
                    gpos1 = gpos2
                else:
                    gpos1 = affected_blocks.end.values[0]
                    gpos2 = gpos1

                novel_seq_info = (chrom, gpos1, gpos2, gene_out.strand.values[0])
                block_seqs['%s:%d-%d(%s)' % novel_seq_info] = str(novel_seq)
                block = pd.DataFrame([[chrom, gpos1, gpos2, 'novel', '.', gene_strand,
                                       'novel_junction', 'novel']], columns=gene_out.columns)
                gene_out = gene_out.append(block, ignore_index=True)
            else:
                block_seqs['%s:%d-%d(%s)' % novel_seq_info] = str(novel_seq)
                block = pd.DataFrame([[chrom, gpos1, gpos2, 'novel', '.', gene_strand,
                                       'unknown', 'novel']], columns=gene_out.columns)
                gene_out = gene_out.append(block, ignore_index=True)
        else:
            block_seqs['%s:%d-%d(%s)' % novel_seq_info] = str(novel_seq)
            blocks_affected = gene_df[np.logical_and(gene_df.start < gpos1, gene_df.end > gpos2)]
            if len(blocks_affected) > 0:
                block = blocks_affected.loc[blocks_affected.index.values[0]]
                gene_out, block_seqs = split_block(nc_row, gene_out, seq, block, gpos1, gpos2)

        gene_out = gene_out.drop_duplicates().sort_values(by=['start', 'end'], ascending=False).reset_index(drop=True) \
                        if antisense else gene_out.drop_duplicates().sort_values(by=['start','end']).reset_index(drop=True)
        # only write normal gene if reference supertranscript is unmodified
        if not gene_out.equals(gene_df):
            var_count += 1
            var = '_v%s' % var_count if len(nc) > 1 else ''
            write_gene(gene_out, block_seqs, '%s_%s%s' % (gene, sample, var), st_file)
            gene_out = gene_df.copy()

    # write modified supertranscript
    write_gene(gene_df, block_seqs, gene, st_file)
