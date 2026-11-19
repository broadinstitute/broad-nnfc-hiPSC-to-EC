#!/usr/bin/env bash
# The v2 guide design, realigned with --tolerate-5p 1.
#
# 45 guides came back no_site from the v2 run because the designer wrote a G
# at position 1 of their protospacer where the genome holds another base, so
# an exact 20-base search found nothing and element-mapping called them
# non-targeting. --tolerate-5p 1 searches the 3'-most 19 bases exactly and
# lets that base differ; the protospacer, its coordinates and the accession
# are unchanged. Nothing else differs from 2026-11-12: the inputs are that
# run's, copied into this one's inputs/, and every other flag is the same.
#
# Run from anywhere; cau comes from the crispr-analysis-utils pixi dev
# environment.
set -euo pipefail

REPO="$(cd "$(dirname "$0")/../.." && pwd)"
CAU_REPO="${CAU_REPO:-$HOME/GitHub/crispr-analysis-utils}"
OUT="$REPO/results/2026-11-18_guide_design_v3"
IN="$OUT/inputs"
STORE="$HOME/Annotations/HG38/GENCODE/gencode.v43.annotation.parquet"
# The same reference as 2026-11-12, md5 a6da8681616c05eb542f1d91606a7b2f.
REFERENCE="${REFERENCE:-$HOME/GitHub/cloud9-mer/reference/hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fa}"

cau() { pixi run --manifest-path "$CAU_REPO/pyproject.toml" -e dev cau "$@"; }

cau guide-alignment run \
    --guides "$REPO/metadata/guides.tsv" \
    --reference "$REFERENCE" \
    --tolerate-5p 1 \
    --outdir "$OUT"

cau guide-annotation run \
    --alignment "$OUT/alignments/guide_genome_annotation.tsv" \
    --store "$STORE" \
    --outdir "$OUT"

cau element-mapping run \
    --alignment "$OUT/alignments/guide_genome_annotation.tsv" \
    --annotation "$OUT/annotations/guide_gene_annotation.tsv" \
    --elements "$IN/elements_combined_v2.bed" \
    --genes "$IN/positive_control_genes_v2.tsv" \
    --outdir "$OUT"

cau igvf-export run \
    --alignment "$OUT/alignments/guide_genome_annotation.tsv" \
    --mapping "$OUT/elements/guide_element_mapping.tsv" \
    --elements "$IN/elements_combined_v2.bed" \
    --store "$STORE" \
    --genes "$IN/positive_control_genes_v2.tsv" \
    --outdir "$OUT"

# cau pairs is not on main yet; it runs from the feat/pairs branch:
#
#     cau pairs run \
#         --guide-metadata "$OUT/igvf/guide_metadata.tsv" \
#         --genes "$REPO/metadata/metadata_genes_final.csv" \
#         --store "$STORE" \
#         --max-distance 1000000 \
#         --outdir "$OUT"
