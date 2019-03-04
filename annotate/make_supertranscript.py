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
import logging
import os
from pybedtools import BedTool
from pybedtools.featurefuncs import extend_fields
from Bio import SeqIO
from argparse import ArgumentParser
from utils import cached, init_logging, exit_with_error
import ipdb

pd.set_option("mode.chained_assignment", None)

EXIT_FILE_IO_ERROR = 1
PROGRAM_NAME = 'MAKE_SUPERTRANSCRIPT'

# headers for GTF file
GTF_COLS = ['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

# only these variant types require modification to reference supertranscripts
VARS_TO_ANNOTATE = ['EE','NE','INS','RI','UN','FUS']

# regex masks
STRAND = '\(([+-])\)'

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
    parser.add_argument(dest='outdir',
                        metavar='OUTDIR',
                        type=str,
                        help='''Output directory.''')
    parser.add_argument(dest='sample',
                        metavar='SAMPLE',
                        type=str,
                        help='''Sample name. Used to name bed and supertranscript output files.''')

    return parser.parse_args()

def reverse_complement(seq):
    # currently unused method, may require STs to be RC'd in the future
    lookup = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
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

def add_strand(exon, strand):
    '''
    add strand to exon block
    '''
    exon = extend_fields(exon, 6)
    exon.strand = strand
    return exon

def get_output_files(sample, outdir):

    genome_bed = '%s/%s_genome.bed' % (outdir, sample)
    st_block_bed = '%s/%s_blocks_supertranscript.bed' % (outdir, sample)
    st_gene_bed = '%s/%s_genes_supertranscript.bed' % (outdir, sample)
    st_fasta = '%s/%s_supertranscript.fasta' % (outdir, sample)

    return genome_bed, st_block_bed, st_gene_bed, st_fasta

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

def split_block(blocks, block, block_seqs, gpos1, gpos2, seq, name, strand):
    '''
    split a sequence block at the gpos location, separating
    into left and right blocks and sequences (if necessary)
    '''
    blocks = blocks[blocks.index!=block['index']]
    ref_seq = block_seqs['%s:%d-%d(%s)' % (block['chr'], block.start, block.end, strand)]
    left_seq = ref_seq[:gpos1-block.start]
    right_seq = ref_seq[gpos2-block.start-1:]

    block_seqs['%s:%d-%d(%s)' % (block['chr'], gpos1, gpos2, strand)] = str(seq)
    var_block = [{'chr': block['chr'], 'start': gpos1, 'end': gpos2, \
                  'name': name, 'score': '.', 'strand': strand}]
    blocks = blocks.append(var_block, ignore_index=True)

    if gpos1 - block.start != 0:
        block_seqs['%s:%d-%d(%s)' % (block['chr'], block.start, gpos1, strand)] = left_seq
        left_block = [{'chr': block['chr'], 'start': block.start, 'end': gpos1, \
                       'name': block['name'], 'score': '.', 'strand': strand}]
        blocks = blocks.append(left_block, ignore_index=True)

    if block.end - gpos2 != 0:
        block_seqs['%s:%d-%d(%s)' % (block['chr'], gpos2, block.end, strand)] = right_seq
        right_block = [{'chr': block['chr'], 'start': gpos2, 'end': block.end, \
                        'name': block['name'], 'score': '.', 'strand': strand}]
        blocks = blocks.append(right_block, ignore_index=True)

    return blocks, block_seqs

def get_merged_exons(genes, gtf, genome_fasta, strand):
    '''
    get all exons from specified genes, merging any
    overlapping exonic regions, also return their
    respective sequences in a dictionary object
    '''
    gene_gtf = gtf[gtf.gene.isin(genes)]
    gene_gtf = gene_gtf.drop('gene', axis=1)
    strand = gene_gtf.strand.values[0] if strand == '' else strand

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
    if strand == '-':
        block_names = ['|' + str(i) for i in reversed(range(1, len(blocks)+1))]
    else:
        block_names = ['|' + str(i) for i in range(1, len(blocks)+1)]
    blocks['name'] = blocks['name'] + block_names

    return(blocks, block_seqs)

def write_supertranscript_genes(blocks, block_bed, gtf, genes, st_gene_bed):
    '''
    write gene annotation to created supertranscript reference
    '''
    seg_starts, seg_ends = block_bed.start, block_bed.end
    contig_name = block_bed['chr'].values[0]

    gene_gtf = gtf[gtf.feature == 'gene']
    if len(gene_gtf) == 0:
        aggregator = {'start': lambda x: min(x),
                      'end': lambda x: max(x)}
        gene_gtf = gtf.groupby(['chr', 'gene', 'strand'], as_index=False, sort=False).agg(aggregator)

    gene_names, gene_starts, gene_ends, gene_strands = [], [], [], []
    for gene in genes:
        gn = gene_gtf[gene_gtf.gene == gene]
        if len(gn) == 0:
            logging.info('WARNING: gene %s not found in reference GTF' % gene)
            continue

        start, end = gn.start.values[0] - 1, gn.end.values[0]
        start_block = blocks[np.logical_and(blocks.start <= start, blocks.end >= start)]
        end_block = blocks[np.logical_and(blocks.start <= end, blocks.end >= end)]
        start_offset = start - min(start_block.start)
        end_offset = max(end_block.end) - end

        gene_start = seg_starts[start_block.index[0]] + start_offset
        gene_end = seg_ends[end_block.index[0]] - end_offset
        gene_starts.append(gene_start)
        gene_ends.append(gene_end)

        # relative to the ST, only assign antisense direction
        # if the reference strand differs from the block strand
        block_strand = blocks[blocks.name.str.contains(gene)].strand.values[0]
        ref_strand = gn.strand.values[0]
        gene_strand = '+' if block_strand == ref_strand else '-'
        gene_strands.append(gene_strand)

        gene_names.append(gene)

    #TODO: add random colours for genes
    if len(gene_starts) > 0:
        bed = pd.DataFrame({'chr': contig_name, 'start': gene_starts, 'end': gene_ends,
                            'name': gene_names, 'score': '.', 'strand': gene_strands})
        bed.to_csv(st_gene_bed, mode='a', index=False, header=False, sep='\t')

def get_block_colours(blocks, names):
    colours = np.empty((len(blocks),), dtype='U50')
    colours[::2] = COL1
    colours[1::2] = COL2
    novel_vars = [x in VARS_TO_ANNOTATE for x in names]
    colours[novel_vars] = VARCOL
    return colours

def write_gene(contig, blocks, block_seqs, args, genes, gtf):
    '''
    write the supertranscript fasta and bed files
    for the given gene. Passing an empty string into
    the contig argument assumes that the gene to write
    is canonical (gene is unmodified from reference)
    '''
    sample = args.sample
    genome_bed, st_block_bed, st_gene_bed, st_fasta = get_output_files(args.sample, args.outdir)

    seqs = []
    for idx,x in blocks.iterrows():
        seq = str(block_seqs['%s:%d-%d(%s)' % (x['chr'], x.start, x.end, x.strand)])
        seqs.append(seq)

    seg_ends = np.cumsum([len(s) for s in seqs])
    seg_starts = np.concatenate([[0], seg_ends[:-1]])
    segs = ['%s-%s' % (s1+1, s2) for s1,s2 in zip(seg_starts, seg_ends)]

    genes = [gn for gn in genes if gn != '']
    names = blocks['name'].apply(lambda x: x.split('|')[-1]).values
    contig_name = '%s|%s|%s' % (sample, contig, '|'.join(genes)) if contig != '' else genes[0]
    header = '>%s segs:%s names:%s\n' % (contig_name, ','.join(segs), ','.join(names))

    # write supertranscript fasta
    sequence = ''.join(seqs) + '\n'
    with open(st_fasta, 'a') as st_fasta:
        st_fasta.writelines([header, sequence])

    # write supertranscript block bed annotation
    colours = get_block_colours(blocks, names)
    bed = pd.DataFrame({'chr': contig_name, 'start': seg_starts, 'end': seg_ends,
                        'name': names, 'score': 0, 'strand': '.', 'thickStart': seg_starts,
                        'thickEnd': seg_ends, 'itemRgb': colours})
    bed.to_csv(st_block_bed, mode='a', index=False, header=False, sep='\t')

    # write supertranscript gene bed annotation
    write_supertranscript_genes(blocks, bed, gtf, genes, st_gene_bed)

def get_contig_genes(con_info):
    '''
    return gene1 and gene2 (in case of fusion),
    indicating genes overlapping the given contig
    '''
    fus_genes = con_info[con_info.overlapping_genes.str.contains(':')]
    if len(fus_genes) > 0:
        fus_genes = np.unique(fus_genes.overlapping_genes)
        fus_genes = [fg.split(':') for fg in fus_genes][0]
        return fus_genes[0], fus_genes[1]
    else:
        genes = np.unique(con_info.overlapping_genes.values)
        if len(genes) > 1:
            logging.info('WARNING: multiple overlapping genes found for contig %s' % con_info.contig_id.values[0])
        return genes[0], ''

def get_contig_strand(con_info, variant):
    '''
    return contig alignment strand given the variant ID
    '''
    strand = '.'
    if variant in con_info.variant_id.values:
        var_info = con_info[con_info.variant_id == variant]
        strand = re.search(STRAND, var_info.pos1.values[0]).group(1)
    if variant in con_info.partner_id.values:
        var_info = con_info[con_info.partner_id == variant]
        strand = re.search(STRAND, var_info.pos2.values[0]).group(1)
    return strand

def sort_blocks(blocks):
    '''
    sort blocks in ascending order if on the sense strand,
    and descending order if on the antisense strand
    '''
    blocks = blocks.drop_duplicates().sort_values(by=['start','end']).reset_index(drop=True)
    if '-' in blocks.strand.values:
        antisense_blocks = blocks[blocks.strand=='-'].sort_values(by=['start','end'], ascending=False)
        sense_blocks = blocks[blocks.strand=='+']
        if len(sense_blocks) > 0:
            sense_first = blocks.strand.values[0] == '+'
            blocks = sense_blocks.append(antisense_blocks) \
                        if sense_first \
                        else antisense_blocks.append(sense_blocks)
        else:
            blocks = antisense_blocks
    return blocks

def contig_to_supertranscript(con_info, args, cvcf, gtf):
    '''
    take the contig info and VCF outputs from annotate/refine contig annotations
    and output three files:
        1) a bed file containing the genomic coordinates of the blocks that
           comprise the gene's supertranscript
        2) a fasta file containing the sequence of each gene's supertranscript
        3) a bed file containing the block annotations (mapping to the ST)
           with novel bits indicated
    '''
    contig = con_info.contig_id.values[0]

    genome_fasta = args.fasta
    genome_bed, st_block_bed, st_gene_bed, st_fasta = get_output_files(args.sample, args.outdir)

    convars = con_info[con_info.variant_type != 'FUS']
    convars = list(convars.variant_id.values) + list(convars.partner_id.values)
    convars = np.unique([c for c in convars if c != '.'])

    is_fusion = 'FUS' in con_info.variant_type.values
    if len(convars) == 0 and not is_fusion:
        # single gene with no variants to annotate
        return

    strands = [re.search(STRAND, con_info.pos1.values[0]).group(1)]
    if 'FUS' in con_info.variant_type.values:
        con_fus = con_info[con_info.variant_type == 'FUS']
        s1 = get_contig_strand(con_fus, con_fus.variant_id.values[0])
        s2 = get_contig_strand(con_fus, con_fus.partner_id.values[0])
        strands = [s1, s2]

    genes = get_contig_genes(con_info)
    blocks, block_seqs = pd.DataFrame(), {}
    for gene, strand in zip(genes, strands):
        if gene != '':
            gene_blocks, gene_block_seqs = get_merged_exons(gene.split('|'), gtf, genome_fasta, strand)
            blocks = blocks.append(gene_blocks)
            for gbs in gene_block_seqs:
                block_seqs[gbs] = gene_block_seqs[gbs]

    if len(blocks) == 0 | len(blocks.drop_duplicates()) != len(block_seqs):
        return

    vcf_records = cvcf[cvcf[2].apply(lambda x: x in convars)] if len(convars) > 0 else pd.DataFrame()
    for idx,record in vcf_records.iterrows():
        chrom = record[0]

        vtype = re.search('SVTYPE=(\w+)', record[7])
        vtype = vtype.group(1) if vtype else 'UN'

        seq = re.search('([ATGCNatgc]+)', record[4])
        if not seq:
            continue
        seq = seq.group(1)
        seq = seq[1:] if vtype == 'INS' else seq
        strand = get_contig_strand(con_info, record[2])
        seq = reverse_complement(seq) if strand == '-' else seq

        blocksize = len(seq) if vtype in ['EE', 'NE', 'RI'] else 0
        start_pos = int(record[1])
        end_pos = int(start_pos) + 1 if blocksize == 0 else int(start_pos) + blocksize

        name = '|'.join(genes) + '|' + vtype
        block_affected = blocks[np.logical_and(blocks.start < start_pos, blocks.end > start_pos)]

        # minor coordinate correction for left-sided soft-clips
        left_sc = vtype == 'UN' and bool(re.search(']', record[4]))
        start_pos = (start_pos - 1) if left_sc else start_pos
        end_pos = (end_pos - 1) if left_sc else end_pos

        if vtype in ['INS', 'UN'] and len(block_affected) > 0:
            if len(block_affected) > 1:
                logging.info('WARNING: multiple blocks affected by variant; exons may not have been merged properly')
            block = block_affected.reset_index().loc[0]
            blocks, block_seqs = split_block(blocks, block, block_seqs, start_pos, end_pos, seq, name, strand)
        else:
            blocks = blocks.append([{'chr': chrom, 'start': start_pos, 'end': end_pos, 'name': name, 'score': '.', 'strand': strand}])
            block_seqs['%s:%d-%d(%s)' % (chrom, start_pos, end_pos, strand)] = seq

    logging.info('Writing contig %s' % contig)
    blocks = sort_blocks(blocks)
    blocks.to_csv(genome_bed, mode='a', index=False, header=False, sep='\t')
    genes = genes[0] + '|' + genes[1] if genes[1] != '' else genes[0]
    write_gene(contig, blocks, block_seqs, args, genes.split('|'), gtf)

def make_supertranscripts(args, contigs, cvcf, gtf):
    '''
    wrapper function for making ST annotation from contigs
    '''
    contigs_to_annotate = contigs[contigs.variant_type.apply(lambda x: x in VARS_TO_ANNOTATE)]
    contig_ids = np.unique(contigs_to_annotate.contig_id.values)

    con_infos = []
    for contig in contig_ids:
        con_info = contigs_to_annotate[contigs_to_annotate.contig_id == contig]
        contig_to_supertranscript(con_info, args, cvcf, gtf)

def write_canonical_genes(args, contigs, gtf):
    '''
    append unmodified reference genes for competitive mapping
    '''
    genes = contigs.overlapping_genes.apply(lambda x: x.split('|'))
    genes = [g.split(':') for gene in genes for g in gene]
    genes = [g for gene in genes for g in gene if g != '']
    genes = np.unique(np.array(genes))

    for idx, gene in enumerate(genes):
        logging.info('Writing %s, %s of %s canonical genes' % (gene, idx+1, len(genes)))
        blocks, block_seqs = get_merged_exons([gene], gtf, args.fasta, '')
        if len(blocks) == 0 | len(blocks) != len(block_seqs):
            continue
        blocks = sort_blocks(blocks)
        write_gene('', blocks, block_seqs, args, [gene], gtf)

def get_gene(attribute):
    '''
    extract gene name from a single GTF attribute string
    '''
    re_gene = re.search('gene_name "([\w\-\.\/]+)"', attribute)
    gene = re_gene.group(1) if re_gene else ''
    return gene

def load_gtf_file(gtf_file):
    '''
    load in reference GTF file containing gene/exon info
    remove 'chr' prefix if present and extract gene names
    '''
    gtf = pd.read_csv(gtf_file, comment='#', sep='\t', header=None, names=GTF_COLS)

    # no non-standard chroms will be handled
    # TODO: is there some way to properly handle alt contigs?
    alt_chrs = gtf['chr'].str.contains('Un|alt|unknown|random|K')
    gtf = gtf[np.invert(alt_chrs.values)]

    # extract gene name from gtf and remove 'chr' prefix if present
    gtf_chrs = gtf['chr'].str.contains('chr')
    if any(gtf_chrs.values):
        gtf['chr'] = gtf['chr'].apply(lambda a: a.split('chr')[1])
        gtf.loc[gtf['chr'] == 'M', 'chr'] = 'MT'
    gtf['gene'] = gtf.attribute.apply(lambda x: get_gene(x))

    return gtf

def load_vcf_file(contig_vcf):
    '''
    load in VCF file containing novel contig variants
    remove 'chr' prefix from chroms if present
    '''
    cvcf = pd.read_csv(contig_vcf, sep='\t', header=None, comment='#')
    vcf_chrs = cvcf[0].str.contains('chr')
    if any(vcf_chrs.values):
        cvcf = cvcf[vcf_chrs]
        cvcf[0] = cvcf[0].apply(lambda a: a.split('chr')[1])
        cvcf.loc[cvcf[0] == 'M', 0] = 'MT'
    return cvcf

def main():
    args = parse_args()
    init_logging(args.log)

    genome_bed, st_block_bed, st_gene_bed, st_fasta = get_output_files(args.sample, args.outdir)
    if os.path.exists(genome_bed):
        os.remove(genome_bed)
    if os.path.exists(st_block_bed):
        os.remove(st_block_bed)
    if os.path.exists(st_gene_bed):
        os.remove(st_gene_bed)
    if os.path.exists(st_fasta):
        os.remove(st_fasta)

    try:
        gtf = load_gtf_file(args.gtf_file)
        cvcf = load_vcf_file(args.contig_vcf)
        contigs = pd.read_csv(args.contig_info, sep='\t').fillna('')
    except IOError as exception:
        exit_with_error(str(exception), EXIT_FILE_IO_ERROR)

    make_supertranscripts(args, contigs, cvcf, gtf)
    write_canonical_genes(args, contigs, gtf)

if __name__ == '__main__':
    main()
