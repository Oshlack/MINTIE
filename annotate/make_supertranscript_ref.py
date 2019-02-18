'''
Module      : make_supertranscript_ref
Description : Make reference supertranscript for cryptic variants
Copyright   : (c) Marek Cmero, Sep 2018
License     : TBD
Maintainer  : MAREK.CMERO@MCRI.EDU.AU
Portability : POSIX
Take VCF, contig_info and reference files, and make a supertranscript
for each contig including any novel bits inserted into the ref sequence
'''

import numpy as np
import pandas as pd
import re
import sys
import tempfile
import os
from pybedtools import BedTool
from Bio import SeqIO
from argparse import ArgumentParser
from utils import cached, init_logging, exit_with_error
import ipdb

pd.set_option("mode.chained_assignment", None)

# for reverse-complementing
lookup = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

EXIT_FILE_IO_ERROR = 1

# headers for GTF file
GTF_COLS = ['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

# TODO: add fusions
# only these variant types require modification to reference supertranscripts
VARS_TO_ANNOTATE = ['EE','NE','INS','RI','UN']

# alternating colours for bed track, and variant colour
COL1 = '189,189,189' #light grey
COL2 = '99,99,99' # dark grey
VARCOL = '255,255,153' # bright yellow

def parse_args():
    '''
    Parse command line arguments.
    Returns Options object with command line argument values as attributes.
    Will exit the program on a command line error.
    '''
    description = 'Make supertranscript reference'
    parser = ArgumentParser(description=description)
    parser.add_argument('--log',
                        metavar='LOG_FILE',
                        type=str,
                        help='record program progress in LOG_FILE')
    parser.add_argument(dest='contig_info',
                        metavar='CONTIG_INFO',
                        type=str,
                        help='''Contig information for novel contigs.''')
    parser.add_argument(dest='contig_vcf',
                        metavar='CONTIG_VCF',
                        type=str,
                        help='''Novel variants in VCF format.''')
    parser.add_argument(dest='gtf_file',
                        metavar='GTF_FILE',
                        type=str,
                        help='''GTF annotation file containing transcript annotations.''')
    parser.add_argument(dest='fasta',
                        metavar='FASTA',
                        type=str,
                        help='''Genome reference in fasta format.''')
    parser.add_argument(dest='sample',
                        metavar='SAMPLE',
                        type=str,
                        help='''Sample name. Used to name bed and supertranscript output files.''')

    return parser.parse_args()

def reverse_complement(seq):
    if seq == '':
        return ''
    if type(seq) == float and math.isnan(seq):
        return ''
    seq = seq[::-1]
    seq = ''.join([lookup[base] for base in list(seq)])
    return(seq)

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

def split_block(blocks, block, block_seqs, gpos1, gpos2, seq, genes, variant):
    gene = '|'.join(genes)

    blocks = blocks[blocks.index!=block['index']]
    ref_seq = block_seqs['%s:%d-%d' % (block['chr'], block.start, block.end)]
    left_seq = ref_seq[:gpos1-block.start]
    right_seq = ref_seq[gpos2-block.start:]

    block_seqs['%s:%d-%d' % (block['chr'], gpos1, gpos2)] = str(seq)
    var_block = [{'chr': block['chr'], 'start': gpos1, 'end': gpos2, 'name': '%s %s' % (gene, variant)}]
    blocks = blocks.append(var_block, ignore_index=True)

    if gpos1 - block.start != 0:
        block_seqs['%s:%d-%d' % (block['chr'], block.start, gpos1)] = left_seq
        left_block = [{'chr': block['chr'], 'start': block.start, 'end': gpos1, 'name': block['name']}]
        blocks = blocks.append(left_block, ignore_index=True)

    if block.end - gpos2 != 0:
        block_seqs['%s:%d-%d' % (block['chr'], gpos2, block.end)] = right_seq
        right_block = [{'chr': block['chr'], 'start': gpos2, 'end': block.end, 'name': block['name']}]
        blocks = blocks.append(right_block, ignore_index=True)

    return blocks, block_seqs

def get_merged_exons(genes, gtf, genome_fasta):
    gene_gtf = pd.DataFrame()
    for gene in genes:
        gene_name = 'gene_name "%s"' % gene
        tmp = gtf[gtf.attribute.str.contains(gene_name)]
        gene_gtf = gene_gtf.append(tmp)

    with tempfile.NamedTemporaryFile(mode='r+') as temp_gtf:
        gene_gtf.to_csv(temp_gtf.name, index=False, header=False, sep='\t')

        # load gene GTF info, extract and merge exons
        g = BedTool(temp_gtf.name)
        exons = BedTool(subset_featuretypes(g, 'exon'))
        exons = exons.remove_invalid().sort().merge()

        exons = exons.sequence(fi=genome_fasta)
        block_seqs = get_block_seqs(exons)

        blocks = pd.DataFrame()
        with tempfile.NamedTemporaryFile(mode='r+') as temp_exons:
            exons.saveas(temp_exons.name)
            blocks = pd.read_csv(temp_exons, header=None, sep='\t', names=['chr', 'start', 'end'])

        if type(genes) == str:
            blocks['name'] = genes
        else:
            blocks['name'] = '|'.join(list(genes))

    blocks['chr'] = blocks['chr'].map(str)
    blocks.start = blocks.start.map(int)
    blocks.end = blocks.end.map(int)
    blocks['name'] = blocks['name'] + ['|' + str(i) for i in range(1, len(blocks)+1)]

    return(blocks, block_seqs)

def write_gene(contig, blocks, block_seqs, st_file, st_bed, genes, sample):
    seqs = []
    for idx,x in blocks.iterrows():
        seq = str(block_seqs['%s:%d-%d' % (x['chr'], x.start, x.end)])
        seqs.append(seq)

    seg_ends = np.cumsum([len(s) for s in seqs])
    seg_starts = np.concatenate([[0], seg_ends[:-1]])
    segs = ['%s-%s' % (s1+1, s2) for s1,s2 in zip(seg_starts, seg_ends)]

    names = blocks['name'].apply(lambda x: x.split('|')[-1]).values
    contig_name = '%s|%s|%s' % (sample, contig, '|'.join(genes))
    header = '>%s segs:%s names:%s\n' % (contig_name, ','.join(segs), ','.join(names))

    sequence = ''.join(seqs) + '\n'
    with open(st_file, 'a') as st_fasta:
        st_fasta.writelines([header, sequence])

    # create and output supertranscript bed file
    colours = np.empty((len(blocks),), dtype='U50')
    colours[::2] = COL1
    colours[1::2] = COL2
    novel_vars = [x in VARS_TO_ANNOTATE for x in names]
    colours[novel_vars] = VARCOL

    bed = pd.DataFrame({'chr': contig_name, 'start': seg_starts, 'end': seg_ends,
                        'name': names, 'score': 0, 'strand': '.', 'thickStart': seg_starts,
                        'thickEnd': seg_ends, 'itemRgb': colours})
    bed.to_csv(st_bed, mode='a', index=False, header=False, sep='\t')

def make_supertranscripts(args, contigs, cvcf, gtf):
    genome_fasta = args.fasta
    # output files
    genome_bed = '%s_genome.bed' % args.sample
    st_bed = '%s_supertranscript.bed' % args.sample
    st_file = '%s_supertranscript.fa' % args.sample

    # remove 'chr' from gtf and vcfs
    # TODO: need to check whether chr is actually present...
    gtf['chr'] = gtf['chr'].apply(lambda a: a.split('chr')[1])
    gtf.loc[gtf['chr'] == 'M', 'chr'] = 'MT'

    cvcf[0] = cvcf[0].apply(lambda a: a.split('chr')[1])
    cvcf.loc[cvcf[0] == 'M', 0] = 'MT'

    contigs_to_annotate = contigs[contigs.variant_type.apply(lambda x: x in VARS_TO_ANNOTATE)]
    contig_ids = np.unique(contigs_to_annotate.contig_id.values)

    for idx,contig in enumerate(contig_ids):
        print('%d of %d contigs' % (idx+1, len(contig_ids)))
        con_info = contigs_to_annotate[contigs_to_annotate.contig_id == contig]

        genes = con_info.overlapping_genes.apply(lambda x: x.split('|'))
        genes = np.unique(np.array([y for x in genes.values for y in x]))

        blocks, block_seqs = get_merged_exons(genes, gtf, genome_fasta)
        if len(blocks) == 0:
            continue

        convars = list(con_info.variant_id.values)
        convars.extend(list(con_info.partner_id.values))
        convars = np.unique([c for c in convars if c != '.'])

        vcf_records = cvcf[cvcf[2].apply(lambda x: x in convars)]
        for idx,record in vcf_records.iterrows():
            chrom = record[0]
            seq = re.search('([ATGCNatgc]+)', record[4]).group(1)[1:]

            vtype = re.search('SVTYPE=(\w+)', record[7])
            vtype = vtype.group(1) if vtype else 'UN'

            blocksize = len(seq) if vtype in ['EE', 'NE', 'RI'] else 0
            start_pos = int(record[1])
            end_pos = int(start_pos) + 1 if blocksize == 0 else int(start_pos) + blocksize

            block_affected = blocks[np.logical_and(blocks.start < start_pos, blocks.end > start_pos)]
            if vtype in ['INS', 'UN'] and len(block_affected) > 0:
                if len(block_affected) > 1:
                    print('WARNING: multiple blocks affected by variant; exons may not have been merged properly')
                block = block_affected.reset_index().loc[0]
                blocks, block_seqs = split_block(blocks, block, block_seqs, start_pos, end_pos, seq, genes, vtype)
            else:
                name = '|'.join(genes) + '|' + vtype
                blocks = blocks.append([{'chr': chrom, 'start': start_pos, 'end': end_pos, 'name': name}])
                block_seqs['%s:%d-%d' % (chrom, start_pos, end_pos)] = seq

        blocks = blocks.drop_duplicates().sort_values(by=['start','end']).reset_index(drop=True)
        blocks.to_csv(genome_bed, mode='a', index=False, header=False, sep='\t')

        write_gene(contig, blocks, block_seqs, st_file, st_bed, genes, args.sample)

def main():
    args = parse_args()
    init_logging(args.log)

    genome_bed = '%s_genome.bed' % args.sample
    st_bed = '%s_supertranscript.bed' % args.sample
    st_file = '%s_supertranscript.fa' % args.sample
    if os.path.exists(genome_bed):
        os.remove(genome_bed)
    if os.path.exists(st_bed):
        os.remove(st_bed)
    if os.path.exists(st_file):
        os.remove(st_file)

    try:
        cvcf = pd.read_csv(args.contig_vcf, sep='\t', header=None, comment='#')
        gtf = pd.read_csv(args.gtf_file, comment='#', sep='\t', header=None, names=GTF_COLS)
        contigs = pd.read_csv(args.contig_info, sep='\t')
    except IOError as exception:
        exit_with_error(str(exception), EXIT_FILE_IO_ERROR)

    make_supertranscripts(args, contigs, cvcf, gtf)

if __name__ == '__main__':
    main()
