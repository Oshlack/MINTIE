# Representing cryptic variants in VCF

The VCF v4.2 specification can be found [here](https://samtools.github.io/hts-specs/VCFv4.2.pdf) and contains ways of representing structural variants that can be adapted for representing transcriptomic/cryptic variants.


The VCF format requires a file format header (containing the VCF version) and a header denoting CHROM, POS, ID, REF and ALT alleles, QUAL and FILTER status. The INFO and FORMAT fields are optional, but are required for representing cryptic variants.

Our version of the VCF looks like this:

```
##fileformat=VCFv4.2
##INFO=<ID=EVENT,Number=1,Type=String,Description="ID of event associated to breakend">
##INFO=<ID=PARID,Number=1,Type=String,Description="ID of partner breakend">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Structural variant type">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT  SAMPLE
chr1    36385897    k49_82102a  TGCGGCGCC   ]TGCGGCGCC  .   .   EVENT=k49_82102;PARID=.;SVLEN=0;SVTYPE=EE   GT  .
chr1    43367998    k49_222045a C   C]chr1:43365632]    .   .   EVENT=k49_222045;PARID=k49_222045b;SVLEN=0;SVTYPE=FUS   GT  .
chr1    43365631    k49_222045b C   C]chr1:43367999]    .   .   EVENT=k49_222045;PARID=k49_222045a;SVLEN=0;SVTYPE=FUS   GT  .
```

### Fixed fields

The first six fields are mostly straight-forward:

1. CHROM -- the genomic chromosome/contig on which the cryptic variant occurs.
2. POS -- start position of the cryptic variant.
3. ID -- the ID of the cryptic variant. This is a unique identifier and is based on the contig name with a latter added to the end. For example, the first variant annotated on contig `k49_11083` will be `k49_11083a`.
4. REF -- this contains the reference sequence.
5. ALT -- this contains the variant sequence. Note that because the genomic sequence can be the same (i.e. splice variants do not cause sequence changes), I add some notation to this field to represent splicing alterations.
6. QUAL -- variant quality. This is currently not implemented, but should be at some point.

### INFO and FORMAT fields

The INFO fields are one-per variant, while the FORMAT field allows a value to be set per sample. Currently, the genotype (GT) is the only sample-specific feature, but this might change. The genotype field gives an allelic presence of the variant, where 0 corresponds to the reference, 1 to the first allele in ALT, 2 to the second allele in ALT and so on. These values are separated by '/' or '|'. For example, a heterozygous cryptic variant would have 0/1 genotype. This functionality hasn't been coded into MINTIE yet, but it should be included at some point.

The INFO fields are described as follows:

1. EVENT -- this is the name of the contig which contains the cryptic variant, this will look something like: `k49_11083`. They could be many variants with the same event (many variants on the same contig).
2. PARID -- this contains the ID of the partner variant. For example, fusions require two records, fusion position 1 and fusion position 2. These will be linked by their partners, e.g. a fusion on `k49_11083` might be represented by variants `k49_11083a` and `k49_11083b` which are at two distant loci. This is consistent with how INFO keys are used for structural variants.
3. SVLEN -- the difference in the ALT and REF allele sequence lengths (might be zero).
4. SVTYPE -- may be:
    1. AS (alternative splicing) -- spliced contig alignment where both junctions are known reference exon boundaries, but spliced reads across these boundaries have not previously been observed.
    2. DEL (deletion) -- aligned section of contig that does not contain part of the reference sequence.
    3. EE (extended exon) -- contig alignment extends past an exon boundary.
    4. FUS (fusion) -- contig involving hard-clipping (one part of the contig aligns to one part of the genome, and the other part aligns elsewhere). Fusions do not have to involved genes, although in our filtering, we require at least one gene to be involved in all cryptic variants.
    5. INS (insertion) -- aligned section of the contig that is not in the reference. This could be a tandem duplication, but is not annotated as such yet. Finer annotation is required to call such events.
    6. NE (novel exon) -- an aligned block that does not cross any existing exons. May be intronic or intergenic.
    7. NEJ (novel exon junction) -- spliced contig where neither junction boundary matches the reference.
    8. PNJ (partial novel junction) -- spliced contig where one junction matches a known boundary (the other side is unknown).
    9. RI (retained intron) -- a contig block that spans a whole intron.
    10. UN (unknown) -- by default, all soft-clipped contigs are classified as unknown at the moment. These variants need to be better annotated.

### Representing variants using the ALT field

Everything is annotated against the _genome_.

#### alternative splicing/novel junctions

These variants would be represented like so (REF and ALT respectively) across two records:

```
G   G[chr21:14064526[
A   ]chr21:14083950]A
```

Where the reference base is the last reference base before the junction if the splicing occurs on the right, otherwise the first reference base is taken for the left side of a junction. The 'G[<loc>[' is the same syntax as for representing SVs, and declares the other end of the junction in reference to the REF base. Technically, we can flip the brackets to represent reverse complements (e.g. `G]chr21:14064526]`), but for now, MINTIE will assume everything occurs on the forward strand as the pipeline is not strand-aware. The second record simply shows the other end, with the reference base following the location.

#### fusions

These variant types are annotated for any contig that is hard-clipped. The way these are annotated is very similar to junction annotation, except that the two record may be 'clipped' on the same side.  For example:

```
C   C]chr1:43365632]
C   C]chr1:43367999]
```

This simple indicates that the loci chr1:43365632 and chr1:43367999 are joined by this variant. If the strands are different across the fusion boundary (as is the case in the example), the brackets will be flipped in the 'reverse complement' position (they will face the reference base).

#### retained intron

Are represented in a single record as follows:

```
CACGTGGGTCGGGGA ]CACGTGGGTCGGGGA[
```

The brackets indicate exon boundaries. The reference sequence will likely be the same as the ALT allele, but it does not have to be (the contig might have a SNP or INDEL for instance). This is the case for all retained intron/extended or novel exon variants.

#### extended exon

Similar to retained intro, extended exons are represented in a single record as follows:

```
GAGGTCACGACGGCGTCCGC    GAGGTCACGACGGCGTCCGC[
```

Where the bracket shows the exon boundary (may be on either side).

#### novel exon

Are represented as follows:

```
CAGGCTGATCACTGTGCCTCCAACTCCGTCGTTCCTGTTCCGATCCCCA   [CAGGCTGATCACTGTGCCTCCAACTCCGTCGTTCCTGTTCCGATCCCCA]
```

Where the brackets indicate the ends of the novel exon.

#### insertions

Are represented as standard in VCF:

```
G   GTGGTGTC
```

#### deletions

Also represented as standard:

```
ACCAATTTGTGCCTGCCCCAGATAG   A
```

#### unknown

These variant types are contigs that are soft-clipped against the reference genome. Further annotation of these is required, but for now, they are represented as:

```
C   CGCCCCCCCCCCCCCCCCCCCCCACCAACC]chr7:5527294]C
```

The sequence on the left is the soft-clip sequence, and the location specifies where the soft-clip starts. A soft-clip may be on the left or right side.

### Contig info file

The output of MINTIE also includes an tab-separated file that specifies contig information for contigs that contain cryptic variants. In principle, this information could be included in the VCF file, but is separated to keep the 'variant' and 'contig' as two separate entities. The tsv file contains the following fields:

1. contig_id: ID of the contig (e.g. `k49_11083`)
2. variant_id: ID of the variant ID described (from VCF)
3. partner_id: partner ID of the variant ID from VCF (if any)
4. pos1: genome position in format \<chr>:\<loc\>
5. pos2: second genome position in format \<chr\>:\<loc\> (may be same as pos1)
6. varsize: variant size (from VCF)
7. cpos: position on the contig where the variant occurs
8. contig_varsize: size of the variant _on the contig_
9. contig_len: total length of the contigs
10. contig_cigar: contig's CIGAR string
11. variant_type: as in VCF file
12. overlapping_genes: genes that overlap the variant position(s), where multiple genes are spannned (or they overlap), the genes are separated with a '|'. Where a fusion (hard-clipped contig) is involved, genes are separated with a ':'
