#!/usr/bin/env python3

import sys
from collections import defaultdict

if len(sys.argv) != 2:
    sys.stderr.write("Usage: python debug_peak_pairs.py peaks_x_guides.tsv\n")
    sys.exit(1)

infile = sys.argv[1]

peak_info = {}
peak_to_guides = defaultdict(set)

with open(infile) as f:
    for line in f:
        if not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")

        peak_chr = fields[0]
        peak_start = int(fields[1])
        peak_end = int(fields[2])
        guide_id = fields[6]

        peak_id = peak_chr + ":" + str(peak_start) + "-" + str(peak_end)

        peak_info[peak_id] = (peak_chr, peak_start, peak_end)
        peak_to_guides[peak_id].add(guide_id)

peaks_by_chr = defaultdict(list)
for peak_id, (chrom, start, end) in peak_info.items():
    peaks_by_chr[chrom].append((start, end, peak_id))

print("\t".join([
    "peak_a",
    "peak_b",
    "n_guides_a",
    "n_guides_b",
    "unique_a",
    "unique_b",
    "decision"
]))

for chrom, peaks in peaks_by_chr.items():
    peaks.sort()

    for i in range(len(peaks)):
        start_i, end_i, peak_a = peaks[i]
        guides_a = peak_to_guides[peak_a]

        j = i + 1
        while j < len(peaks) and peaks[j][0] < end_i:
            start_j, end_j, peak_b = peaks[j]
            guides_b = peak_to_guides[peak_b]

            unique_a = len(guides_a - guides_b)
            unique_b = len(guides_b - guides_a)

            if unique_a < 5 and unique_b < 5:
                decision = "merge"
            else:
                decision = "separate"

            print("\t".join(map(str, [
                peak_a,
                peak_b,
                len(guides_a),
                len(guides_b),
                unique_a,
                unique_b,
                decision
            ])))

            j += 1