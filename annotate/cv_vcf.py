import string

CIGAR = {'match': 0,
         'insertion': 1,
         'deletion': 2,
         'skipped': 3,
         'soft-clip': 4,
         'hard-clip': 5,
         'silent_deletion': 6}
AFFECT_CONTIG = [CIGAR['insertion'], CIGAR['match'], CIGAR['soft-clip'], CIGAR['hard-clip']]

# VCF parameters
VCF_VERSION = "4.2"
COLUMNS = ["#CHROM",  "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]

# number, type and description for header output
info_fields = {"SVTYPE": (1, "String", "Structural variant type"),
               "SVLEN": (1, "Integer", "Difference in length between REF and ALT alleles"),
               "PARID": (1, "String", "ID of partner breakend"),
               "EVENT": (1, "String", "ID of event associated to breakend")}
format_fields = {"GT": (1, "String", "Genotype")}

# contig info parameters
CI_HEADER = ["contig_id", "variant_id", "partner_id", "pos1", "pos2",
             "varsize", "cpos", "contig_varsize", "contig_len", "contig_cigar",
             "variant_type", "overlapping_genes"]

class VCF(object):
    '''
    VCF object
    '''
    def __init__(self):
        pass

    @staticmethod
    def get_header(sample):
        print('##fileformat=VCFv%s' % VCF_VERSION)

        # INFO lines
        for info_id in sorted(info_fields):
            info = info_fields[info_id]
            line = '##INFO=<ID=%s,Number=%s,Type=%s,Description="%s">' % (info_id, info[0], info[1], info[2])
            print(line)

        # FORMAT lines
        for format_id in sorted(format_fields):
            form = format_fields[format_id]
            line = '##FORMAT=<ID=%s,Number=%s,Type=%s,Description="%s">' % (format_id, form[0], form[1], form[2])
            print(line)

        # columns line
        line = COLUMNS + [sample]
        return "\t".join(line)

    @staticmethod
    def format_record(chrom, pos, vid, ref, alt, qual, filterval, info, form):
        info_out = ['%s=%s' % (key, str(info[key])) for key in sorted(info_fields)]
        info_out = ';'.join(info_out)

        form_fields = ':'.join(sorted(format_fields))
        form_out = ':'.join([str(form[key]) for key in sorted(format_fields)])

        outline = [chrom, pos, vid, ref, alt, qual, filterval, info_out, form_fields, form_out]
        outline = [str(item) for item in outline]

        return "\t".join(outline)

class CrypticVariant(object):
    '''
    Cryptic contig object
    '''
    def __init__(self, **kwargs):
        "Build an empty CrypticVariant object"
        defaults = {
            "chrom": 'NA',
            "pos": 0,
            "strand": '.',
            "vid": '.',
            "ref": '.',
            "alt": '.',
            "qual": '.',
            "cfilter": '.',
            "cid": '.',
            "clen": 0,
            "cpos": 0,
            "cstrand": '.',
            "ccigar": '.',
            "vsize": 0,
            "cvsize": 0,
            "cvtype": '.',
            "genes": '.',
            "parid": '.',
            "ecn": ".",
            "ai": [0.00, 0.00, 0.00],
            "pval": [1.0, 1.0],
            "cvq": 0.00,
            "gt": ".",
            "ecc": 0,
            "blocks": []
        }
        for (prop, default) in defaults.items():
            setattr(self, prop, kwargs.get(prop, default))

    @staticmethod
    def write_contig_header(contig_info_file):
        '''
        Write header for contig info output
        '''
        with open(contig_info_file, 'w') as fout:
            fout.write("\t".join(CI_HEADER) + "\n")

    @staticmethod
    def write_contig_info(contig_info_file, cv1, cv2=None):
        '''
        Write contig info for a given variant/variant pair
        '''
        # make varsize block size for block annotations
        varsize = cv1.cvsize - 1 if cv1.cvtype in ['NE', 'EE', 'RI'] else cv1.vsize

        pos1 = "%s:%d(%s)" % (cv1.chrom, cv1.pos, cv1.cstrand)
        pos2 = "%s:%d(%s)" % (cv2.chrom, cv2.pos, cv2.cstrand) if cv2 else "%s:%d(%s)" % (cv1.chrom, (cv1.pos + varsize), cv1.cstrand)

        genes = '%s:%s' % (cv1.genes, cv2.genes) if cv2 else cv1.genes
        # only consider a gene its own variant partner if it's a fusion
        if cv1.cvtype != 'FUS' and cv2 and cv1.genes == cv2.genes:
            genes = cv1.genes

        line = [cv1.cid, cv1.vid, cv1.parid, pos1,
                pos2, cv1.vsize, cv1.cpos, cv1.cvsize,
                cv1.clen, cv1.ccigar, cv1.cvtype, genes]
        line = [str(item) for item in line]

        with open(contig_info_file, 'a') as fout:
            fout.write("\t".join(line) + "\n")

    def from_read(self, read):
        match_idxs = [idx for idx,cig in enumerate(read.cigar) if cig[0] == 0]
        self.chrom = read.reference_name
        self.vid = read.query_name + 'a'
        self.cstrand = '-' if read.is_reverse else '+'
        self.blocks = [b for b in zip(match_idxs, read.get_blocks())]
        self.cid = read.query_name
        self.ccigar = read.cigarstring
        self.clen =  sum([v for c,v in read.cigar if c in AFFECT_CONTIG])
        return self

    def get_format(self):
        form = {"GT": self.gt}
        return form

    def get_info(self):
        info = {"SVTYPE": self.cvtype,
                "SVLEN": self.vsize,
                "PARID": self.parid,
                "EVENT": self.cid}
        return info

    def vcf_output(self):
        '''
        Take cryptic variant properties and
        output a VCF formatted record.
        '''
        info, form = self.get_info(), self.get_format()
        return VCF.format_record(self.chrom,
                                 self.pos,
                                 self.vid,
                                 self.ref,
                                 self.alt,
                                 self.qual,
                                 self.cfilter,
                                 info, form)
