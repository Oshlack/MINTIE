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
import logging
import os
import bedtool_helper
import tempfile
import pickle
import block_helper as bh
from pybedtools import BedTool
from Bio import SeqIO
from argparse import ArgumentParser
from utils import cached, init_logging, exit_with_error

pd.set_option("mode.chained_assignment", None)

EXIT_FILE_IO_ERROR = 1

# headers for GTF file
GTF_COLS = ['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
BED_EXT_COLS = ['chr', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb']

# only these variant types require modification to reference supertranscripts
VARS_TO_ANNOTATE = ['EE','NE','INS','RI','UN','FUS','DEL']

# regex masks
STRAND = r'\(([+-])\)'

# keep track of canonical genes written to avoid duplicate entries
canonical_genes_written = []

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

#=====================================================================================================
# Utility functions
#=====================================================================================================

def reverse_complement(seq):
    lookup = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    if seq == '':
        return ''
    if type(seq) == float and math.isnan(seq):
        return ''
    seq = seq[::-1]
    seq = ''.join([lookup[base.upper()] for base in list(seq)])
    return(seq)

def get_gene(attribute):
    '''
    extract gene name from a single GTF attribute string
    '''
    try:
        return attribute['gene_name']
    except KeyError:
        return ''

def get_contig_genes(con_info):
    '''
    return gene1 and gene2 (in case of fusion),
    indicating genes overlapping the given contig
    '''
    fus_genes = con_info[con_info.overlapping_genes.str.contains(':')]
    if len(fus_genes) > 0:
        fus_genes = np.unique(fus_genes.overlapping_genes)
        fus_genes = [fg.split(':') for fg in fus_genes][0]
        assert len(fus_genes) <= 2
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

def get_gene_strands(gtf, genes):
    strands = []
    for gene in genes:
        gene = gene.split('|')[0] # take first gene as representative in case of overlaps
        strand = gtf[gtf.gene == gene].strand.values
        strand = strand[0] if gene != '' and len(strand) > 0 else ''
        strands.append(strand)
    return strands

def get_strand_info(con_info, gstrands):
    if 'FUS' in con_info.variant_type.values:
        con_fus = con_info[con_info.variant_type == 'FUS']
        cs1 = get_contig_strand(con_fus, con_fus.variant_id.values[0])
        cs2 = get_contig_strand(con_fus, con_fus.partner_id.values[0])

        # if one fragment aligns on a strand counter to gene orientation
        # while the other aligns matching gene orientation, return strands
        # corresponding to contig alignment
        gs1, gs2 = gstrands
        if (cs1 != gs1 and cs2 == gs2) or (cs1 == gs1 and cs2 != gs2):
            return [cs1, cs2]
        else:
            return gstrands
    else:
        return gstrands

#=====================================================================================================
# Read/write functions
#=====================================================================================================

def get_output_files(sample, outdir):
    genome_bed = '%s/%s_genome.bed' % (outdir, sample)
    st_block_bed = '%s/%s_blocks_supertranscript.bed' % (outdir, sample)
    st_gene_bed = '%s/%s_genes_supertranscript.bed' % (outdir, sample)
    st_fasta = '%s/%s_supertranscript.fasta' % (outdir, sample)

    return genome_bed, st_block_bed, st_gene_bed, st_fasta

@cached('gene_gtf.pickle')
def load_gtf_file(gtf_file):
    '''
    load in reference GTF file containing gene/exon info
    remove 'chr' prefix if present and extract gene names
    '''
    logging.info('Processing GTF reference file...')

    gtf = BedTool(gtf_file).remove_invalid().saveas()
    gene_names = [get_gene(row.attrs) for row in gtf]

    with tempfile.NamedTemporaryFile(mode='r+') as temp_gtf:
        gtf.saveas(temp_gtf.name)
        gtf_pd = pd.read_csv(temp_gtf, header=None, sep='\t', names=GTF_COLS, comment='#')
        gtf_pd['gene'] = gene_names

    # no non-standard chroms will be handled
    # TODO: is there some way to properly handle alt contigs?
    alt_chrs = gtf_pd['chr'].str.contains('Un|alt|unknown|random|K')
    gtf_pd = gtf_pd[np.invert(alt_chrs.values)]

    # extract gene name from gtf_pd and remove 'chr' prefix if present
    gtf_pd_chrs = gtf_pd['chr'].str.contains('chr')
    if any(gtf_pd_chrs.values):
        gtf_pd['chr'] = gtf_pd['chr'].apply(lambda a: a.split('chr')[1])
        gtf_pd.loc[gtf_pd['chr'] == 'M', 'chr'] = 'MT'

    # create gene features if none exist
    if len(gtf_pd[gtf_pd.feature == 'gene']) == 0:
        aggregator = {'start': lambda x: min(x),
                      'end': lambda x: max(x)}
        gene_gtf_pd = gtf_pd.groupby(['chr', 'score', 'strand', 'frame', 'gene'],
                                     as_index=False, sort=False).agg(aggregator)
        gene_gtf_pd['source'] = 'ALL'
        gene_gtf_pd['attribute'] = ''
        gene_gtf_pd['feature'] = 'gene'
        gene_gtf_pd = gene_gtf_pd[GTF_COLS + ['gene']]
        gtf_pd = gtf_pd.append(gene_gtf_pd)

    return gtf_pd

def load_vcf_file(contig_vcf):
    '''
    load in VCF file containing novel contig variants
    remove 'chr' prefix from chroms if present
    '''
    cvcf = pd.read_csv(contig_vcf, sep='\t', header=None, comment='#')
    vcf_chrs = cvcf[0].map(str).str.contains('chr')
    if any(vcf_chrs.values):
        cvcf = cvcf[vcf_chrs]
        cvcf[0] = cvcf[0].apply(lambda a: a.split('chr')[1])
        cvcf.loc[cvcf[0] == 'M', 0] = 'MT'
    return cvcf

def write_supertranscript_genes(blocks, block_bed, gtf, genes, st_gene_bed):
    '''
    Write gene locations on ST in bed file
    '''
    seg_starts, seg_ends = block_bed.start, block_bed.end
    seg_starts.index, seg_ends.index = blocks.index, blocks.index
    contig_name = block_bed['chr'].values[0]

    gene_gtf = gtf[gtf.feature == 'gene']
    gene_names, gene_starts, gene_ends, gene_strands, gene_cols = [], [], [], [], []

    # get gene strands from block annotation
    # (inferring strand if gene not present)
    block_strands = []
    last_strand = blocks.strand.values[0]
    for gene in genes:
        for gn in gene:
            gn_match = blocks.name.str.contains(gn)
            if any(gn_match):
                last_strand = blocks[gn_match].strand.values[0]
                block_strands.append(last_strand)
            else:
                block_strands.append(last_strand)
    genes = [g for gn in genes for g in gn]

    for gene, block_strand in zip(genes, block_strands):
        gn = gene_gtf[gene_gtf.gene == gene]
        if len(gn) == 0:
            logging.info('WARNING: gene %s not found in reference GTF.' % gene)
            continue

        start, end = gn.start.values[0] - 1, gn.end.values[0]
        start_block = blocks[np.logical_and(blocks.start <= start, blocks.end >= start)]
        end_block = blocks[np.logical_and(blocks.start <= end, blocks.end >= end)]
        if len(start_block) == 0 or len(end_block) == 0:
            logging.info('WARNING: part of gene %s sits outside the reference blocks; skipping.')
            continue
        start_offset = start - min(start_block.start)
        end_offset = max(end_block.end) - end

        # relative to the ST, only assign antisense direction
        # if the reference strand differs from the block strand
        ref_strand = gn.strand.values[0]
        gene_strand = '+' if block_strand == ref_strand else '-'
        gene_strands.append(gene_strand)

        antisense = block_strand == '-'
        tmp = start_block.copy()
        start_block = end_block if antisense else start_block
        end_block = tmp if antisense else end_block
        gene_start = seg_starts[start_block.index[0]] + start_offset
        gene_end = seg_ends[end_block.index[0]] - end_offset

        R, G, B = np.random.randint(0, 255, 3)
        gene_col = '%s,%s,%s' % (R, G, B)

        gene_starts.append(gene_start)
        gene_ends.append(gene_end)
        gene_names.append(gene)
        gene_cols.append(gene_col)

    if len(genes) > 1 and genes[0] == genes[1]:
        gene_starts[1] = gene_ends[0]
        gene_ends[1] = gene_ends[1] * 2

    if len(gene_starts) > 0:
        bed = pd.DataFrame({'chr': contig_name, 'start': gene_starts, 'end': gene_ends,
                            'name': gene_names, 'score': '.', 'strand': gene_strands,
                            'thickStart': gene_starts, 'thickEnd': gene_ends,
                            'itemRgb': gene_cols},
                            columns=BED_EXT_COLS)
        bed.to_csv(st_gene_bed, mode='a', index=False, header=False, sep='\t')

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

    # for contig name, extract first gene as representative
    gene_name = [gn.split('|')[0] for gn in genes if gn != '']
    contig_name = '%s|%s|%s' % (sample, contig, '|'.join(gene_name)) if contig != '' else gene_name[0]
    names = blocks['name'].apply(lambda x: x.split('|')[-1]).values
    header = '>%s segs:%s names:%s\n' % (contig_name, ','.join(segs), ','.join(names))

    # write supertranscript fasta
    sequence = ''.join(seqs) + '\n'
    with open(st_fasta, 'a') as st_fasta:
        st_fasta.writelines([header, sequence])

    # alternate block colour tracks in cases of fusion
    if len(genes) > 1 and '' not in genes:
        blocks = blocks.reset_index()

        gnames = blocks.name.apply(lambda x: x.split('|')[0])
        gns, idxs = np.unique(gnames.values, return_index=True)

        if len(gns) == 1:
            half_idx = int((len(blocks)/2))
            g1 = blocks.index < half_idx
            g2 = blocks.index >= half_idx
        else:
            # get unique genes in the order they appear
            gene1, gene2 = gns[:2] if idxs[0] < idxs[1] else gns[1::-1]
            g1 = gnames == gene1
            g2 = gnames == gene2

        c1 = bh.get_block_colours(blocks[g1], names[g1])
        c2 = bh.get_block_colours(blocks[g2], names[g2], alt=True)
        colours = np.concatenate([c1, c2])

        assert len(colours) == len(seg_starts)
    else:
        colours = bh.get_block_colours(blocks, names)

    # write supertranscript block bed annotation
    bed = pd.DataFrame({'chr': contig_name, 'start': seg_starts, 'end': seg_ends,
                        'name': names, 'score': 0, 'strand': '.', 'thickStart': seg_starts,
                        'thickEnd': seg_ends, 'itemRgb': colours}, columns=BED_EXT_COLS)
    bed.to_csv(st_block_bed, mode='a', index=False, header=False, sep='\t')

    # write supertranscript gene bed annotation
    genes = [gn.split('|') for gn in genes if gn != '']
    write_supertranscript_genes(blocks, bed, gtf, genes, st_gene_bed)

def write_canonical_genes(args, contigs, gtf):
    '''
    append unmodified reference genes for competitive mapping
    '''
    genes = contigs.overlapping_genes.apply(lambda x: x.split('|'))
    genes = [g.split(':') for gene in genes for g in gene]
    genes = [g for gene in genes for g in gene if g != '' and g not in canonical_genes_written]
    genes = np.unique(np.array(genes))

    logging.info('%d additional canonical genes to write...' % len(genes))
    for gene in genes:
        logging.info('Writing %s' % gene)
        blocks, block_seqs = bedtool_helper.get_merged_exons([gene], gtf, args.fasta, '')
        if len(blocks) == 0:
            continue
        if len(blocks.drop_duplicates()) != len(block_seqs):
            continue
        blocks = bh.sort_blocks(blocks)
        write_gene('', blocks, block_seqs, args, [gene], gtf)

#=====================================================================================================
# Processing functions
#=====================================================================================================

def get_block_info(args, genes, strands, gtf, genome_fasta):
    '''
    get blocks and corresponding sequence info
    for all genes overlapping given contig
    '''
    blocks, block_seqs = pd.DataFrame(), {}
    for gene, strand in zip(genes, strands):
        if gene != '':
            gene_blocks, gene_block_seqs = bedtool_helper.get_merged_exons(gene.split('|'),
                                                                           gtf, genome_fasta, strand)
            if len(gene_blocks) == 0:
                 return blocks, block_seqs
            if len(gene_blocks.drop_duplicates()) != len(gene_block_seqs):
                 return blocks, block_seqs
            blocks = blocks.append(gene_blocks)
            for gbs in gene_block_seqs:
                block_seqs[gbs] = gene_block_seqs[gbs]
            if gene not in canonical_genes_written and gene.find('|') < 0:
                logging.info('Writing %s' % gene)
                canonical_genes_written.append(gene)
                write_gene('', bh.sort_blocks(blocks), block_seqs, args, [gene], gtf)
    return blocks.drop_duplicates(), block_seqs

def add_novel_sequence(blocks, block_seqs, record, con_info, genes, strand):
    '''
    add any novel sequence from the given record to the reference blocks
    '''
    chrom = record[0]
    vtype = re.search(r'SVTYPE=(\w+)', record[7])
    vtype = vtype.group(1) if vtype else 'UN'

    seq = re.findall('([ATGCNatgc]+)', record[4]) # TODO: make this parsing more robust
    if not seq:
        return blocks, block_seqs
    seq = max(seq, key=len)
    seq = seq[1:] if vtype in ['INS', 'UN', 'DEL'] else seq
    seq = reverse_complement(seq) if strand == '-' else seq

    blocksize = len(seq) if vtype in ['EE', 'NE', 'RI'] else 0
    start_pos = int(record[1])
    end_pos = int(start_pos) + 1 if blocksize == 0 else int(start_pos) + blocksize
    end_pos = int(start_pos) + (len(record[3]) - len(seq)) if vtype == 'DEL' \
                                                           else end_pos

    name = genes[0] + '|' + vtype
    block_affected = blocks[np.logical_and(blocks.start < start_pos, blocks.end > start_pos)]

    # minor coordinate correction for left-sided soft-clips
    left_sc = vtype == 'UN' and bool(re.search(']', record[4]))
    start_pos = (start_pos - 1) if left_sc else start_pos
    end_pos = (end_pos - 1) if left_sc else end_pos

    if vtype in ['INS', 'UN', 'DEL'] and len(block_affected) > 0:
        if len(block_affected) > 1:
            logging.info('''WARNING: multiple blocks affected by variant;
                            Exons may not have been merged properly''')
        block = block_affected.reset_index().loc[0]
        blocks, block_seqs = bh.split_block(blocks, block, block_seqs, start_pos,
                                            end_pos, seq, name, strand)
    else:
        blocks = blocks.append([{'chr': chrom,
                                 'start': start_pos,
                                 'end': end_pos,
                                 'name': name,
                                 'score': '.',
                                 'strand': strand}])
        block_seqs['%s:%d-%d(%s)' % (chrom, start_pos, end_pos, strand)] = seq
    return blocks, block_seqs

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
    genome_fasta = args.fasta
    genome_bed, st_block_bed, st_gene_bed, st_fasta = get_output_files(args.sample, args.outdir)

    convars = con_info[con_info.variant_type != 'FUS']
    convars = list(convars.variant_id.values) + list(convars.partner_id.values)
    convars = np.unique([c for c in convars if c != '.'])

    is_fusion = 'FUS' in con_info.variant_type.values
    if len(convars) == 0 and not is_fusion:
        # single gene with no variants to annotate
        return

    genes = get_contig_genes(con_info)
    gene_strands = get_gene_strands(gtf, genes)
    strands = get_strand_info(con_info, gene_strands)

    blocks, block_seqs = get_block_info(args, genes, strands, gtf, genome_fasta)
    if len(blocks) == 0 or len(blocks) != len(block_seqs):
        return

    vcf_records = cvcf[cvcf[2].isin(convars)] if len(convars) > 0 else pd.DataFrame()
    for idx,record in vcf_records.iterrows():
        gene = con_info[con_info.variant_id == record[2]].overlapping_genes.values[0]
        strand = [s for g,s in zip(genes, strands) if g == gene][0]
        blocks, block_seqs = add_novel_sequence(blocks, block_seqs, record, con_info, genes, strand)

    contig = con_info.contig_id.values[0]
    logging.info('Writing supertranscript for contig %s' % contig)
    blocks = bh.sort_blocks(blocks)
    blocks.to_csv(genome_bed, mode='a', index=False, header=False, sep='\t')

    if genes[0].split('|')[0] == genes[1].split('|')[0]:
        #TODO: handle case where genes need to be orientated differently or
        #cases where genes overlap in genes 1 and 2 but have a different genes
        #in their respective leading positions
        blocks = blocks.append(blocks, ignore_index=True)

    write_gene(contig, blocks, block_seqs, args, genes, gtf)

def make_supertranscripts(args, contigs, cvcf, gtf):
    '''
    wrapper function for making ST annotation from contigs
    '''
    annotate_contigs = np.logical_and(contigs.variant_of_interest,
                                      contigs.variant_type.isin(VARS_TO_ANNOTATE))
    contigs_to_annotate = contigs[annotate_contigs]
    contig_ids = np.unique(contigs_to_annotate.contig_id.values)

    con_infos = []
    for contig in contig_ids:
        con_info = contigs_to_annotate[contigs_to_annotate.contig_id == contig]
        contig_to_supertranscript(con_info, args, cvcf, gtf)

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
