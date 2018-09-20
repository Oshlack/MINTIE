import string
import re

# CUT-OFF parameters
GAP_MIN = 7
CLIP_MIN = 30
MATCH_MIN = 30
MATCH_PERC_MIN = 0.3

# VCF parameters
INFO = ["CID", "ECN", "CLEN", "CPOS", "CSTRAND", "CCIGAR", "VSIZE", "CVSIZE", "CVTYPE", "GENES", "PARID", "PVAL", "CVQ"]
FORMAT = ["GT", "ECC", "AI"]

# CIGAR specification codes
CIGAR = {'match': 0,
         'insertion': 1,
         'deletion': 2,
         'skipped': 3,
         'soft-clip': 4,
         'hard-clip': 5,
         'silent_deletion': 6}
GAPS = {1: 'insertion',
        2: 'deletion',
        6: 'silent_deletion'}
CLIPS = {4: 'soft',
         5: 'hard'}
GAPS_REF_ONLY = {2: 'deletion',
                 3: 'skipped',
                 5: 'hard-clip'} # criteria that signify gaps in the reference

#TODO: add all metadata and make header generation function
def get_next_letter(last_letter):
    next_letter_pos = np.where(np.array(list(string.ascii_letters)) == last_letter)[0][0]+1
    next_letter = list(string.ascii_letters)[next_letter_pos]
    return next_letter

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
            "cpos": [0,0],
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

    def from_read(self, read):
        match_idxs = [idx for idx,cig in enumerate(read.cigar) if cig[0] == 0]
        self.chrom = read.reference_name
        self.vid = read.query_name + 'a'
        self.cstrand = '-' if read.is_reverse else '+'
        self.blocks = [b for b in zip(match_idxs, read.get_blocks())]
        self.cid = read.query_name
        self.ccigar = read.cigarstring
        self.clen =  sum([v for c,v in read.cigar if c in [0, 1, 4, 5]]) # count if match, insertion or clip TODO: parametrise
        return self

    def get_format(self):
        gt = str(self.gt)
        ecc = str(self.ecc)
        ai = ','.join([str(a) for a in self.ai])
        return ':'.join([gt, ecc, ai])

    def get_info(self):
        cid = str(self.cid)
        ecn = str(self.ecn)
        clen = str(self.clen)
        cpos = ','.join([str(cp) for cp in self.cpos])
        cstrand = str(self.cstrand)
        ccigar = str(self.ccigar)
        vsize = str(self.vsize)
        cvsize = str(self.cvsize)
        cvtype = str(self.cvtype)
        genes = str(self.genes)
        parid = str(self.parid)
        pval = ','.join([str(pv) for pv in self.pval])
        cvq = str(self.cvq)
        output = zip(INFO, [cid, ecn, clen, cpos, cstrand, ccigar, vsize, cvsize, cvtype, genes, parid, pval, cvq])
        return ';'.join(["%s=%s" % (c, v) for c,v in output])

    def vcf_output(self):
        chrom = str(self.chrom)
        pos = str(self.pos)
        vid = str(self.vid)
        ref = str(self.ref)
        alt = str(self.alt)
        qual = str(self.qual)
        cfilter = str(self.cfilter)
        info = self.get_info()
        form = self.get_format()
        return "\t".join([chrom, pos, vid, ref, alt, qual, cfilter, info, ':'.join(FORMAT), form])
