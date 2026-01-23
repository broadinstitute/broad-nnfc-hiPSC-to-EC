import re
import sys
import pandas as pd

# Read your TSV
df = pd.read_csv(sys.argv[1], sep="\t")

# Function to compute reference length from CIGAR
def cigar_ref_len(cigar):
    # Count ops that consume reference: M, D, N, =, X
    length = 0
    for length_str, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):
        if op in "MDN=X":
            length += int(length_str)
    return length

# Compute start (0-based) and end
df["start"] = df["pos1"] - 1
df["end"] = df["start"] + df["cigar"].apply(cigar_ref_len)

# Build BED dataframe
bed_df = df[[
    "chromosome",  # chrom
    "start",
    "end",
    "read_name",   # name
    "mapq",        # score
    "strand",      # strand
    "reason"       # extra column, BED allows extra fields
]]

# Write BED
bed_df.to_csv("output.bed", sep="\t", header=False, index=False)

