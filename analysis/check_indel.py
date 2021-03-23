"""
Exploring the reads supporting a given indel from a BAM file
"""
import datetime
import os
import sys

import pysam

# Input indel
INDEL_CONTEXT = sys.argv[1]
INDEL = ':'.join(INDEL_CONTEXT.split(':')[:-1])
ANNOTATION = sys.argv[2]

# Out files
OUT_PREFIX = INDEL.replace(':', '_')
OUT_DIR = OUT_PREFIX
os.makedirs(OUT_DIR, exist_ok=True)
OUT_SUMMARY_FILE_PATH = os.path.join(
    OUT_DIR, f"{OUT_PREFIX}_summary.txt"
)
OUT_SUMMARY_FILE = open(OUT_SUMMARY_FILE_PATH, 'w')
OUT_SUMMARY_FILE.write(f"{datetime.date.today().isoformat()}\n")
OUT_REF_ALLELE_READS_FILE_PATH = os.path.join(
    OUT_DIR, f"{OUT_PREFIX}_ref_allele.txt"
)
OUT_REF_ALLELE_READS_FILE = open(OUT_REF_ALLELE_READS_FILE_PATH, 'w')
OUT_ALT_ALLELE_READS_FILE_PATH = os.path.join(
    OUT_DIR, f"{OUT_PREFIX}_alt_allele.txt"
)
OUT_ALT_ALLELE_READS_FILE = open(OUT_ALT_ALLELE_READS_FILE_PATH, 'w')

# General statistics
OUT_SUMMARY_FILE.write(f"Explored indel: {INDEL}\t{ANNOTATION}\n")
CHR = INDEL.split(':')[0]
START = int(INDEL.split(':')[1]) # Position of the indel in 1-base
REF = INDEL.split(':')[2]
ALT = INDEL.split(':')[3]
END = START + len(REF) - 1
# Sequence starting at POS and including the full reference allele
REF_CONTEXT = INDEL_CONTEXT.split(':')[-1]
REF_CONTEXT_LEN = len(REF_CONTEXT)
END_REF_CONTEXT = START + REF_CONTEXT_LEN - 1
OUT_SUMMARY_FILE.write(f"Reference allele context: {REF_CONTEXT}\n")
# Sequence starting at POS and including the full alternate allele
ALT_CONTEXT = f"{ALT}{REF_CONTEXT[len(REF):]}"
ALT_CONTEXT_LEN = len(ALT_CONTEXT)
END_ALT_CONTEXT = START + ALT_CONTEXT_LEN - 1
OUT_SUMMARY_FILE.write(f"Alternate allele context: {ALT_CONTEXT}\n")

# Path to input BAM file
BAM_FILE_PATH = sys.argv[3]
BAM_FILE = pysam.AlignmentFile(BAM_FILE_PATH, 'rb')
OUT_SUMMARY_FILE.write(f"Input BAM file: {BAM_FILE_PATH}\n")

# Auxiliary functions

def read_sign(mapping):
    if mapping.is_read1:
        return '1'
    if mapping.is_read2:
        return '2'

def stats_for_mappings(mappings, description):
    MAPPINGS_NB = len(mappings)
    OUT_SUMMARY_FILE.write(
        f"\nThe {description} contains {MAPPINGS_NB} mappings\n"
    )
    READ_PAIRS_NAMES = set([mapping.query_name for mapping in mappings])
    READ_PAIRS_NB = len(READ_PAIRS_NAMES)
    OUT_SUMMARY_FILE.write(
        f"The {description} originate "
        f"from {READ_PAIRS_NB} distinct read pairs\n"
    )
    READ_NAMES = set(
        [f"{mapping.query_name}:{read_sign(mapping)}" for mapping in mappings]
    )
    READS_NB = len(READ_NAMES)
    OUT_SUMMARY_FILE.write(
        f"The {description} originate from {READS_NB} distinct reads\n"
    )
    return (MAPPINGS_NB, READ_PAIRS_NB, READS_NB)

def check_context(mapping, context_seq):
    query_context_start = START - (mapping.reference_start + 1)
    query_context_end = START - (mapping.reference_start + 1) + len(context_seq)
    query_seq = mapping.query_alignment_sequence
    query_context_seq = query_seq[query_context_start:query_context_end]
    return query_context_seq == context_seq

# Analysis

# Fetching mappings covering the indel
MAPPINGS = list(BAM_FILE.fetch(contig=CHR, start=START, stop=END + 1))
MAPPINGS_NB = stats_for_mappings(
    MAPPINGS, 'full set of overlapping mappings'
)

# Primary mappings
PRIMARY_MAPPINGS = [mapping for mapping in MAPPINGS if not mapping.is_secondary]
stats_for_mappings(PRIMARY_MAPPINGS, 'set of primary overlapping mappings')

# Mappings covering the full reference allele
FULL_MAPPINGS = [
    mapping for mapping in MAPPINGS
    if (mapping.reference_start + 1) <= START and (mapping.reference_end >= END)
]
PRIMARY_MAPPINGS_NB = stats_for_mappings(
    FULL_MAPPINGS, 'set of mappings covering the reference allele'
)

# Mappings covering the full reference allele context
FULL_REF_CONTEXT_MAPPINGS = [
    mapping for mapping in MAPPINGS
    if (
        (mapping.reference_start + 1) <= START and
        (mapping.reference_end >= END_REF_CONTEXT)
    )
]
FULL_REF_CONTEXT_MAPPINGS_NB = stats_for_mappings(
    FULL_REF_CONTEXT_MAPPINGS,
    'set of mappings covering the reference allele context'
)

# Mappings with the reference allele sequence
REF_CONTEXT_MAPPINGS = [
    mapping for mapping in FULL_REF_CONTEXT_MAPPINGS
    if check_context(mapping, REF_CONTEXT)
]
REF_CONTEXT_MAPPINGS_NB = stats_for_mappings(
    REF_CONTEXT_MAPPINGS,
    'set of mappings including the reference allele context'
)
REF_CONTEXT_READS = list(set([
    f"{mapping.query_name}:{read_sign(mapping)}"
    for mapping in REF_CONTEXT_MAPPINGS
]))
REF_CONTEXT_READS.sort()
OUT_REF_ALLELE_READS_FILE.write('\n'.join(REF_CONTEXT_READS))

# Mappings with the alternate  allele sequence
ALT_CONTEXT_MAPPINGS = [
    mapping for mapping in FULL_REF_CONTEXT_MAPPINGS
    if check_context(mapping, ALT_CONTEXT)
]
ALT_CONTEXT_MAPPINGS_NB = stats_for_mappings(
    ALT_CONTEXT_MAPPINGS,
    'set of mappings including the alternate allele context'
)
ALT_ALLELE_READS_FILE = open(f"{INDEL.replace(':', '_')}_alt_allele.txt", 'w')
ALT_CONTEXT_READS = list(set([
    f"{mapping.query_name}:{read_sign(mapping)}"
    for mapping in ALT_CONTEXT_MAPPINGS
]))
ALT_CONTEXT_READS.sort()
OUT_ALT_ALLELE_READS_FILE.write('\n'.join(ALT_CONTEXT_READS))

OUT_REF_ALLELE_READS_FILE.close()
OUT_ALT_ALLELE_READS_FILE.close()
OUT_SUMMARY_FILE.close()
