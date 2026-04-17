#!/usr/bin/env python3

import sys
from collections import defaultdict

if len(sys.argv) != 3:
    sys.stderr.write("Usage: python merge_by_guides_overlap_only.py peaks_x_guides.tsv out.bed\n")
    sys.exit(1)

infile = sys.argv[1]
outfile = sys.argv[2]

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

class UF:
    def __init__(self, items):
        self.parent = {x: x for x in items}

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra = self.find(a)
        rb = self.find(b)
        if ra != rb:
            self.parent[rb] = ra

all_peak_ids = list(peak_info.keys())
uf = UF(all_peak_ids)

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
                uf.union(peak_a, peak_b)

            j += 1

groups = defaultdict(list)
for peak_id in all_peak_ids:
    groups[uf.find(peak_id)].append(peak_id)

rows = []
for members in groups.values():
    members = sorted(
        set(members),
        key=lambda p: (peak_info[p][0], peak_info[p][1], peak_info[p][2], p)
    )

    chrom = peak_info[members[0]][0]
    start = min(peak_info[p][1] for p in members)
    end = max(peak_info[p][2] for p in members)
    guides = sorted(set().union(*(peak_to_guides[p] for p in members)))

    rows.append((
        chrom,
        start,
        end,
        chrom + ":" + str(start) + "-" + str(end),
        len(members),
        len(guides),
        ",".join(members),
        ",".join(guides)
    ))

rows.sort(key=lambda x: (x[0], x[1], x[2]))

with open(outfile, "w") as out:
    for row in rows:
        out.write("\t".join(map(str, row)) + "\n")