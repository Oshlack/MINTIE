'''
Module      : refine_annotations
Description : Perform further filtering on annotated contigs
Copyright   : (c) Marek Cmero, Feb 2019
License     : MIT
Maintainer  : github.com/mcmero
Portability : POSIX
'''

import numpy as np
import pandas as pd
import re
import sys
import logging
import pysam
import annotate_contigs as ac
import pybedtools as pbt
import constants
import tempfile
from Bio import SeqIO
from intervaltree import IntervalTree
from pybedtools import BedTool
from argparse import ArgumentParser
from utils import init_logging, exit_with_error

PROGRAM_NAME = 'refine_annotations'

SPLICE_VARS = ['AS']
SV_VARS = ['DEL', 'INS']
NOVEL_BLOCKS = ['EE', 'NE']
NOVEL_JUNCS = ['PNJ', 'NEJ']
LARGE_SV = ['FUS', 'IGR']
UNKNOWN = ['UN']
SENSE_MOTIF = ['AG', 'GT']
ANTISENSE_MOTIF = ['AC', 'CT']

def parse_args():
    '''
    Parse command line arguments.
    Returns Options object with command line argument values as attributes.
    Will exit the program on a command line error.
    '''
    description = 'Filter contigs'
    parser = ArgumentParser(description=description)
    parser.add_argument('--log',
                        metavar='LOG_FILE',
                        type=str,
                        help='record program progress in LOG_FILE')
    parser.add_argument(dest='contig_info_file',
                        metavar='CONTIG_INFO_FILE',
                        type=str,
                        help='''Contig info file''')
    parser.add_argument(dest='vcf_file',
                        metavar='VCF_FILE',
                        type=str,
                        help='''Contig VCF file''')
    parser.add_argument(dest='bam_file',
                        metavar='BAM_FILE',
                        type=str,
                        help='''Contig BAM file''')
    parser.add_argument(dest='tx_ref_file',
                        type=str,
                        metavar='TX_REF_FILE',
                        help='''Transcriptiome GTF reference file.''')
    parser.add_argument(dest='fasta',
                        metavar='FASTA',
                        type=str,
                        help='''Genome fasta file''')
    parser.add_argument(dest='out_prefix',
                        metavar='OUT_PREFIX',
                        type=str,
                        help='''Output prefix''')
    parser.add_argument('--minClip',
                        metavar='MIN_CLIP',
                        type=int,
                        help='''Minimum novel block or softclip size.''')
    parser.add_argument('--minGap',
                        metavar='MIN_GAP',
                        type=int,
                        help='''Minimum gap (deletion or insertion) size.''')
    parser.add_argument('--mismatches',
                        metavar='MISMATCHES',
                        type=int,
                        default=0,
                        help='''Number of allowed mismatches when checking splice motifs (0-4).
                              0 = no mismatches allowed (default).
                              1 = allows a single mismatch (total) across both the donor and acceptor sites.
                              2 = allows two mismatches (one in both donor and acceptor).
                              3 = motifs are not checked but are returned for certain variants.
                              4 = motifs are not checked or returned.
                              Note that options 1 and 2 work identically when checking a single site only.''')
    return parser.parse_args()

def set_globals(args):
    global MIN_CLIP
    global MIN_GAP

    if args.minClip:
        MIN_CLIP = args.minClip
    else:
        MIN_CLIP = constants.DEFAULT_MIN_CLIP

    if args.minGap:
        MIN_GAP = args.minGap
    else:
        MIN_GAP = constants.DEFAULT_MIN_GAP

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

def load_vcf_file(contig_vcf):
    '''
    load in VCF file containing novel contig variants
    remove 'chr' prefix from chroms if present
    '''
    cvcf = pd.read_csv(contig_vcf, sep='\t', header=None, comment='#', low_memory=False)
    return cvcf

def get_diff_count(motif, side = 0, sense = True):
    '''
    Returns the number of differences (with respect to position)
    from the splicing motifs, given a 2-character motif, side
    and strand.
    '''
    assert len(motif) == 2
    test_motif = SENSE_MOTIF[side] if sense else ANTISENSE_MOTIF[side]
    count = int(motif[0] != test_motif[0])
    count += int(motif[1] != test_motif[1])
    return count

def is_motif_valid(motifs, mismatches):
    '''
    Given a motif (left and right ends), returns whether
    it is valid on either strand within an allowed number
    of mismatches.
    '''
    if mismatches > 2:
        return True
    elif '' in motifs:
        # single-end motif
        side = np.where('' != np.array(motifs))[0][0]
        diff_s = get_diff_count(motifs[side], side = side)
        diff_as = get_diff_count(motifs[side], side = side, sense = False)

        # mismatches == 2 only affects both ends
        mismatches = 1 if mismatches == 2 else mismatches
        return diff_s <= mismatches or diff_as <= mismatches
    else:
        # check motif on both ends
        ldiff_s = get_diff_count(motifs[0])
        ldiff_as = get_diff_count(motifs[0], sense = False)
        rdiff_s = get_diff_count(motifs[1], side = 1)
        rdiff_as = get_diff_count(motifs[1], side = 1, sense = False)

        # check if sense strand valid
        if ldiff_s + rdiff_s <= mismatches and \
            ldiff_s < 2 and rdiff_s < 2:
            return True

        # check if antisense strand valid
        if ldiff_as + rdiff_as <= mismatches and \
            ldiff_as < 2 and rdiff_as < 2:
            return True
    return False

def check_valid_motif(left_id, right_id, block_seqs, mismatches):
    '''
    Test whether motif is valid for given sequences.
    '''
    try:
        lseq = '' if left_id == '' else block_seqs[left_id]
        rseq = '' if right_id == '' else block_seqs[right_id]

        motif = [lseq, rseq]
        valid = False if lseq == '' and rseq == '' else is_motif_valid(motif, mismatches)
        return valid, ''.join(motif)
    except KeyError:
        # occurs if chrom not in reference
        logging.info('''WARNING: one of the following motif Locations could not be
                        retrieved from the provided reference: %s, %s''' % (left_id, right_id))
        return False, ''

def get_valid_motif_vars(variants, args):
    '''
    Checks whether the given variant has a canonical splice site motif,
    returning those variant IDs with valid motifs. Novel exons are checked
    for valid donor and acceptor sites (in either transcriptional direction),
    while extended exons on the left and right are checked for either valid
    donor or acceptor sites (in either transcriptional direction).
    '''
    #TODO: refactor function; too long
    # get VCF data of given variants
    vcf = load_vcf_file(args.vcf_file)
    vcf = vcf[vcf[2].apply(lambda x: x in variants.variant_id.values)]
    vcf = vcf[vcf[3].apply(lambda x: np.invert(pd.isnull(x)))]

    # construct motif locations for which to extract sequence
    b_blocks = vcf[vcf[4].apply(lambda x: x.startswith('[') and x.endswith(']'))]
    r_blocks = vcf[vcf[4].apply(lambda x: x.startswith(']'))]
    l_blocks = vcf[vcf[4].apply(lambda x: x.endswith('['))]
    left = pd.concat([l_blocks, b_blocks])
    right = pd.concat([r_blocks, b_blocks])
    mlocs = pd.DataFrame({'chr': pd.concat([left[0], right[0]]),
                          'start': pd.concat([left[1] - 3, right[1] + right[3].apply(len)-1])})
    mlocs['end'] = mlocs['start'] + 2
    mlocs = mlocs.drop_duplicates()

    # ensure coords aren't negative or exceed chrom len
    mlocs = mlocs[mlocs['start'] >= 0]
    chr_sizes = pbt.chromsizes('hg38') #TODO: allow reference to be argument
    for chrom in np.unique(mlocs['chr'].values):
        ref_chrom = 'chr%s' % chrom if chrom != 'MT' else 'chrM'
        try:
            chr_max = chr_sizes[ref_chrom][1]
            over_limit = np.logical_or(mlocs.start > chr_max, mlocs.end > chr_max)
            mlocs = mlocs.drop(mlocs[np.logical_and(mlocs['chr'] == chrom, over_limit)].index)
        except KeyError:
            logging.info("WARNING: unable to check length for chrom %s, as it does not exist in hg38 reference." % ref_chrom)
            continue

    # extract sequences
    g = BedTool.from_dataframe(mlocs).remove_invalid()
    g = g.sequence(fi=args.fasta)
    bs = get_block_seqs(g)

    # check whether variants have valid motifs
    # left blocks
    valid_vars, motifs = [], []
    if len(l_blocks) > 0:
        left_ids = ['%s:%d-%d' % loc for loc in zip(l_blocks[0], l_blocks[1]-3, l_blocks[1]-1)]
        info_left = [check_valid_motif(lid, '', bs, args.mismatches) for lid in left_ids]
        valid_left = [v for (v, m) in info_left]
        motifs.extend([m for (v, m) in info_left])
        if any(valid_left):
            valid_vars.extend(l_blocks[valid_left][2].values)

    # right blocks
    if len(r_blocks) > 0:
        rpos = r_blocks[1] + r_blocks[3].apply(len)-1
        right_ids = ['%s:%d-%d' % loc for loc in zip(r_blocks[0], rpos, rpos+2)]
        info_right = [check_valid_motif('', rid, bs, args.mismatches) for rid in right_ids]
        valid_right = [v for (v, m) in info_right]
        motifs.extend([m for (v, m) in info_right])
        if any(valid_right):
            valid_vars.extend(r_blocks[valid_right][2].values)

    # both blocks (novel exons)
    if len(b_blocks) > 0:
        rpos = b_blocks[1] + b_blocks[3].apply(len)-1
        left_ids = ['%s:%d-%d' % loc for loc in zip(b_blocks[0], b_blocks[1]-3, b_blocks[1]-1)]
        right_ids = ['%s:%d-%d' % loc for loc in zip(b_blocks[0], rpos, rpos+2)]
        info_both = [check_valid_motif(lid, rid, bs, args.mismatches) for lid, rid in zip(left_ids, right_ids)]
        valid_both = [v for (v, m) in info_both]
        motifs.extend([m for (v, m) in info_both])
        if any(valid_both):
            valid_vars.extend(b_blocks[valid_both][2].values)

    motif_info = pd.concat([l_blocks, r_blocks, b_blocks])
    motif_info['motif'] = motifs
    motif_info['variant_id'] = motif_info[2]
    motif_info['valid_motif'] = motif_info.variant_id.isin(valid_vars)
    return motif_info[['variant_id', 'motif', 'valid_motif']]

def check_for_valid_motifs(contigs, vars_to_check, args):
    if any(vars_to_check):
        motif_info = get_valid_motif_vars(contigs[vars_to_check], args)
        contigs = contigs.merge(motif_info, on = 'variant_id', how = 'left')
        contigs['motif'] = contigs.motif.fillna('')
    return contigs

def check_overlap(ex_trees, chrom, start, end, size=0):
    '''
    Checks whether variant overlaps an exonic region.
    For deletions, at least MIN_GAP bp of the deletion
    must be within the exon body.
    '''
    olap = False
    ex_tree = ac.get_chrom_ref_tree(chrom, ex_trees)
    if not ex_tree:
        return olap

    olap = ex_tree.overlaps(start, end)
    if olap and size > 0:
        olap_se = ex_tree.overlap(start, end)
        es, ee = [(x[0], x[1]) for x in olap_se][0]
        size_within = min([ee, end]) - start if start >= es \
                                      else end - max([es, start])
        olap = size_within >= size

    return olap

def get_pos_parts(loc):
    loc = loc.split(':')
    chrom = loc[0]
    pos = int(loc[1].split('(')[0])
    strand = '.'
    try:
        strand = re.search(r'\(([-+])\)', loc[1]).group(1)
    except AttributeError:
        logging.info('WARNING: invalid strand in %s.' % loc)

    return chrom, pos, strand

def get_varsize(sv):
    '''
    Variant end - start size
    '''
    chr1, start, s1 = get_pos_parts(sv['pos1'])
    chr2, end, s2 = get_pos_parts(sv['pos2'])

    return end - start

def overlaps_same_exon(sv, ex_trees):
    '''
    Checks whether variant is contained
    completely within a single exon
    '''
    chr1, start, s1 = get_pos_parts(sv['pos1'])
    chr2, end, s2 = get_pos_parts(sv['pos2'])

    ex_tree = ac.get_chrom_ref_tree(chr1, ex_trees)
    if ex_tree:
        olap1 = ex_tree.overlap(start, start+1)
        olap2 = ex_tree.overlap(end, end+1)
        return len(olap1) > 0 and len(olap2) > 0 and olap1 == olap2

    return False

def overlaps_exon(sv, ex_trees):
    '''
    Checks whether variant overlaps an exonic region.
    In variants involving two ends (fusions, junctions etc.)
    only one end has to overlap and exon to return true
    '''
    chr1, start, s1 = get_pos_parts(sv['pos1'])
    chr2, end, s2 = get_pos_parts(sv['pos2'])

    # check starts/ends separately for all vars but deletions and juncs
    span_vars = ['DEL'] + NOVEL_JUNCS
    end = end if sv['variant_type'] in span_vars else start + 1

    size = MIN_GAP if sv['variant_type'] == 'DEL' else 0
    size = MIN_CLIP if sv['variant_type'] in NOVEL_JUNCS else size

    olap = check_overlap(ex_trees, chr1, start, end, size=size)
    if sv['variant_type'] in LARGE_SV:
        olap = olap or check_overlap(ex_trees, chr2, end,
                                     start + 1, size=size)
    return olap

def match_splice_juncs(contigs):
    '''
    Check whether exons contain matching novel/partial novel junctions
    '''
    spliced_exons = []
    exons = contigs[contigs.variant_type.isin(NOVEL_BLOCKS)]
    # some junction gaps may be called as deletions, so include these
    novel_juncs = contigs[contigs.variant_type.isin(NOVEL_JUNCS + ['DEL'])]
    for idx,row in exons.iterrows():
        back_junc = novel_juncs.pos2 == row['pos1']
        front_junc = novel_juncs.pos1 == row['pos2']
        matching_juncs = novel_juncs[np.logical_or(back_junc, front_junc)]
        if len(matching_juncs) > 0:
            spliced_exons.append(row['variant_id'])

    is_spliced_exon = contigs.variant_id.isin(spliced_exons)
    return is_spliced_exon

def vars_overlap_exon(contigs, ex_trees):
    '''
    Checks whether variants overlap an exonic region
    '''
    vars_overlaps_exon = np.array([False] * len(contigs))
    exonic_var = contigs.variant_type.isin(['EE', 'AS'])
    contigs.loc[exonic_var, 'overlaps_exon'] = True
    check_overlap = np.invert(exonic_var)
    vars_overlaps_exon[check_overlap] = contigs[check_overlap].apply(overlaps_exon,
                                                                     axis=1, args=(ex_trees,))
    return vars_overlaps_exon

def get_junc_vars(contigs, ex_trees, args):
    '''
    Return truncated exons and novel introns
    '''
    # check for novel exon juncs contained within single exon (may be deletions)
    within_exon = contigs.apply(overlaps_same_exon, axis=1, args=(ex_trees,))
    nj_var = contigs.variant_type.isin(NOVEL_JUNCS)
    nj_dels = np.empty(0, dtype=object)
    if sum(nj_var.values) > 0:
        bigger_than_mingap = contigs[nj_var].apply(get_varsize, axis=1) >= MIN_GAP
        nj_dels = contigs[nj_var][np.logical_and(within_exon[nj_var], bigger_than_mingap)].variant_id.values

    # check truncated-exon vars
    is_trunc = np.logical_and.reduce((contigs.variant_type.isin(NOVEL_JUNCS),
                                      np.invert(within_exon),
                                      contigs.overlaps_exon))
    if 'valid_motif' in contigs.columns.values:
        contigs_tmp = contigs.drop(['motif', 'valid_motif'], axis = 1)
        contigs_tmp = check_for_valid_motifs(contigs_tmp, is_trunc, args)
        contigs.loc[is_trunc, 'motif'] = contigs_tmp.motif
        contigs.loc[is_trunc, 'valid_motif'] = contigs_tmp.valid_motif
        is_trunc = np.logical_and(is_trunc, contigs.valid_motif)
    trunc_vars = contigs[is_trunc].variant_id.values

    return np.unique(np.concatenate([nj_dels, trunc_vars]))

def get_tsv_vars(contigs):
    '''
    TSV criteria:
    - larger than MIN_GAP
    - overlaps exon
    - is softclipped (UN) and larger than MIN_CLIP
    '''
    # check whether TSVs meets size requirements
    is_sv = contigs.variant_type.isin(SV_VARS)
    large_gap = np.logical_or(contigs.varsize >= MIN_GAP,
                              contigs.contig_varsize >= MIN_GAP)
    keep_sv = np.logical_and(large_gap, is_sv)
    keep_sv = np.logical_and(keep_sv, contigs.overlaps_exon)
    sv_vars = contigs[keep_sv].variant_id.values

    # keep unknown (soft-clipped) variants
    is_un = contigs.variant_type.isin(UNKNOWN)
    large_clip = contigs.varsize >= MIN_CLIP
    un_vars = contigs[np.logical_and(is_un, large_clip)].variant_id.values
    sv_vars = np.concatenate([un_vars, sv_vars])

    return sv_vars

def get_fusion_vars(contigs):
    '''
    Return all variants and associated variants at
    large rearrangement (fusion/IGR) boundaries
    '''
    is_fus = contigs.variant_type.isin(LARGE_SV)
    fus_ids = contigs[is_fus].contig_id.values
    fus_locs = np.union1d(contigs[is_fus].pos1, contigs[is_fus].pos2)
    non_fus_vars = contigs[np.logical_and(contigs.contig_id.isin(fus_ids), np.invert(is_fus))]

    # consider variants at boundaries of large rearrangements interesting
    at_fus_boundary = np.logical_or(non_fus_vars.pos1.isin(fus_locs),
                                    non_fus_vars.pos2.isin(fus_locs))
    fus_boundary_vars = non_fus_vars[at_fus_boundary].variant_id.values
    fus_vars = contigs[is_fus].variant_id.values
    fus_vars = np.concatenate([fus_vars, fus_boundary_vars])

    return fus_vars

def overlaps_gene(row, gene_tree):
    '''
    Return True if variant in contig row
    overlaps any gene in the reference
    '''
    olaps = False
    chr1, pos1, strand1 = get_pos_parts(row['pos1'])
    chr2, pos2, strand2 = get_pos_parts(row['pos2'])
    if chr1 == chr2:
        gtree = ac.get_chrom_ref_tree(chr1, gene_tree)
        olaps = gtree.overlaps(pos1, pos2) if gtree else False
    else:
        gtree = ac.get_chrom_ref_tree(chr1, gene_tree)
        olaps = gtree.overlaps(pos1, pos1 + 1) if gtree else False
        gtree = ac.get_chrom_ref_tree(chr2, gene_tree)
        olaps = olaps or gtree.overlaps(pos2, pos2 + 1) if gtree else olaps
    return olaps

def get_contigs_to_keep(args):
    '''
    Return contigs matching criteria:
    - novel block size >= MIN_CLIP
    - novel blocks are spliced in some way
    - novel exons have corresponding novel splice
      sites and novel splice donor/acceptor motifs
    - TSVs occur in exonic regions
    '''
    try:
        cinfo_file = args.contig_info_file
        contigs = pd.read_csv(cinfo_file, sep='\t')
    except IOError as exception:
        exit_with_error(str(exception), constants.EXIT_FILE_IO_ERROR)

    gene_tree, ex_trees, ex_ref = ac.get_gene_lookup(args.tx_ref_file)
    contigs['large_varsize'] = contigs.contig_varsize >= MIN_CLIP
    contigs['is_contig_spliced'] = contigs.contig_cigar.apply(lambda x: bool(re.search('N', x)))
    contigs['spliced_exon'] = match_splice_juncs(contigs)
    contigs['overlaps_exon'] = vars_overlap_exon(contigs, ex_trees)
    contigs['overlaps_gene'] = contigs.apply(overlaps_gene, axis=1, args=(gene_tree,))

    # keep exons falling outside the gene
    is_intergenic_exon = np.logical_and.reduce((contigs.spliced_exon,
                                                contigs.variant_type == 'NE',
                                                contigs.large_varsize,
                                                np.invert(contigs.overlaps_gene)))
    # novel exon contigs (spliced, valid motif and large variant size)
    is_novel_exon = np.logical_and(contigs.spliced_exon, contigs.large_varsize)
    if args.mismatches < 4:
        contigs = check_for_valid_motifs(contigs, is_novel_exon, args)
        is_novel_exon = np.logical_and(is_novel_exon, contigs.valid_motif)
    ne_vars = contigs[np.logical_or(is_novel_exon, is_intergenic_exon)].variant_id.values

    # keep all splice vars
    as_vars = contigs.variant_id.values[contigs.variant_type.isin(SPLICE_VARS)]

    # get junc vars
    junc_vars = get_junc_vars(contigs, ex_trees, args)

    # check size of retained introns
    ri_vars = contigs[np.logical_and(contigs.variant_type == 'RI',
                                     contigs.large_varsize)].variant_id.values

    # get fusions and TSVs
    fus_vars = get_fusion_vars(contigs)
    sv_vars = get_tsv_vars(contigs)

    # collate contigs to keep
    keep_vars = np.unique(np.concatenate([ri_vars, as_vars, ne_vars, sv_vars, fus_vars, junc_vars]))
    contigs['variant_of_interest'] = contigs.variant_id.isin(keep_vars)
    contigs.to_csv('%s_info.tsv' % args.out_prefix, sep='\t', index=None)

    keep_contigs = contigs[contigs.variant_of_interest].contig_id.values
    return(keep_contigs)

def write_output(args, keep_contigs):
    cvars_file = args.vcf_file
    try:
        vcf = pd.read_csv(cvars_file, sep='\t', header=None, comment='#')
        cvf = open(cvars_file, 'r')
        for line in cvf:
            if not line.startswith('#'):
                break
            print(line.strip())
        cvf.close()
    except IOError as exception:
        exit_with_error(str(exception), constants.EXIT_FILE_IO_ERROR)

    vcf = vcf[vcf[7].apply(lambda x: x.split(';')[0].split('=')[1] in keep_contigs)]
    vcf.to_csv(sys.stdout, sep='\t', index=False, header=False)

def write_bam(args, keep_contigs):
    bam_file = args.bam_file
    bam = pysam.AlignmentFile(bam_file, 'rb')
    outbam = pysam.AlignmentFile('%s.bam'% args.out_prefix, 'wb', template=bam)

    for read in bam.fetch():
        if read.query_name in keep_contigs:
            outbam.write(read)

def main():
    args = parse_args()
    init_logging(args.log)
    set_globals(args)
    keep_contigs = get_contigs_to_keep(args)
    if len(keep_contigs) > 0:
        write_output(args, keep_contigs)
        write_bam(args, keep_contigs)
    else:
        exit_with_error('ERROR: no variants to output.', constants.EXIT_OUTPUT_ERROR)

if __name__ == '__main__':
    main()
