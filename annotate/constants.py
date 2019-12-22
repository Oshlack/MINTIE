#=====================================================================================================
# Program parameters
#=====================================================================================================
EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
EXIT_OUTPUT_ERROR = 3

#=====================================================================================================
# Default values for cutoff parameters
#=====================================================================================================
DEFAULT_MIN_GAP = 3
DEFAULT_MIN_CLIP = 30
DEFAULT_MIN_MATCH_BP = 30
DEFAULT_MIN_MATCH_PERC = 0.3

#=====================================================================================================
# VCF output parameters
#=====================================================================================================
INFO = ["CID", "ECN", "CLEN", "CPOS", "CSTRAND", "CCIGAR", "VSIZE",
        "CVSIZE", "CVTYPE", "GENES", "PARID", "PVAL", "CVQ"]
FORMAT = ["GT", "ECC", "AI"]

#=====================================================================================================
# CIGAR-string related
#=====================================================================================================
CIGAR = {'match': 0,
         'insertion': 1,
         'deletion': 2,
         'skipped': 3,
         'soft-clip': 4,
         'hard-clip': 5,
         'silent_deletion': 6}

GAPS = [CIGAR[c] for c in ['insertion', 'deletion', 'silent_deletion']]
CLIPS = [CIGAR[c] for c in ['soft-clip', 'hard-clip']]

# any cigar criteria that is >0 bp on an aligned contig
AFFECT_CONTIG = [CIGAR[c] for c in ['match', 'insertion', 'soft-clip', 'hard-clip']]

# any cigar criteria that is >0 bp on the reference genome
AFFECT_REF = [CIGAR[c] for c in ['match', 'deletion']]
