#!/usr/bin/env python3
import csv

A_BED = "../../results/03-2025-10-24_create_sceptre_input/elements_bed_final_4col.bed"
A_FIELDS = 4          # your A has 4 columns (from your example)
B_FIELDS = 6          # our gencode bed is BED6
KEY_COLS_A = (0, 1, 2, 3)  # chr,start,end,element_id

B_NAME_COL = A_FIELDS + 3          # B "name" is the 4th column of B (0-based +3)
OVERLAP_COL = A_FIELDS + B_FIELDS  # overlap is appended after A+B

def make_keys():
    keys = set()
    with open(A_BED) as f:
        for row in csv.reader(f, delimiter="\t"):
            if row:
                keys.add((row[0], row[1], row[2], row[3]))
    return keys

def read_best_hit(path):
    best = {}  # key -> (overlap, gene_id|gene_name)
    with open(path) as f:
        for row in csv.reader(f, delimiter="\t"):
            if not row:
                continue
            key = tuple(row[i] for i in KEY_COLS_A)
            ov = int(row[OVERLAP_COL])
            gene = row[B_NAME_COL]
            cur = best.get(key)
            if cur is None or ov > cur[0]:
                best[key] = (ov, gene)
    return best

keys = make_keys()

gene_best = read_best_hit("../../results/03-2025-10-24_create_sceptre_input/el_vs_gene.wo.tsv")
exon_best = read_best_hit("../../results/03-2025-10-24_create_sceptre_input/el_vs_exon.wo.tsv")
cds_best  = read_best_hit("../../results/03-2025-10-24_create_sceptre_input/el_vs_cds.wo.tsv")
utr_best  = read_best_hit("../../results/03-2025-10-24_create_sceptre_input/el_vs_utr.wo.tsv")

def pick_feature(key):
    if key in cds_best:  return ("CDS",  cds_best[key][1])
    if key in utr_best:  return ("UTR",  utr_best[key][1])
    if key in exon_best: return ("exon", exon_best[key][1])
    if key in gene_best: return ("intron", gene_best[key][1])
    return ("intergenic", ".")

with open("../../results/03-2025-10-24_create_sceptre_input/elements.annotated.tsv", "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["chr","start","end","element_id","feature","gene_id|gene_name"])
    for key in sorted(keys):
        feature, gene = pick_feature(key)
        w.writerow([key[0], key[1], key[2], key[3], feature, gene])

print("Wrote elements.annotated.tsv")

