"""
Extracting Natera indels that overlap our amplicon panel
"""

import io
import pandas as pd
import sys

INDEL_VCF_INPUT_FILE = sys.argv[1]

def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})
INDELS_INPUT_DF = read_vcf(INDEL_VCF_INPUT_FILE)

MANIFEST_FILE_NAME = 'CG001v5.1_Amplicon_Manifest_Panel5.1.14_20201119.tsv'
MANIFEST_DF = pd.read_csv(MANIFEST_FILE_NAME, sep='\t')
MANIFEST_DF = MANIFEST_DF.astype({'Start': int, 'End': int})

INDELS_OUTPUT_INDEX = []
for _, amplicon in MANIFEST_DF.iterrows():
    start, end = amplicon['Start'], amplicon['End']
    indels_index = INDELS_INPUT_DF.loc[
        (INDELS_INPUT_DF['POS'].between(start, end)) &
        ((len(INDELS_INPUT_DF['REF']) > 1) | (len(INDELS_INPUT_DF['ALT']) > 1))
    ].index
    INDELS_OUTPUT_INDEX += list(indels_index)

VCF_HEADER = (
    f"##fileformat=VCFv4.2"
    f"\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    f"\n##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">"
    f"\n#"
)

INDELS_OUTPUT_DF = INDELS_INPUT_DF.iloc[list(set(INDELS_OUTPUT_INDEX))]

INDELS_VCF_OUTPUT_FILE_NAME = INDEL_VCF_INPUT_FILE.replace(
    '.vcf',
    f"_{MANIFEST_FILE_NAME.replace('.tsv', '.vcf')}"
)
INDELS_VCF_OUTPUT_FILE = open(INDELS_VCF_OUTPUT_FILE_NAME, 'w')
INDELS_VCF_OUTPUT_FILE.write(VCF_HEADER)
INDELS_VCF_OUTPUT_FILE.close()
INDELS_OUTPUT_DF.to_csv(INDELS_VCF_OUTPUT_FILE_NAME, sep='\t', mode='a', index=False)
