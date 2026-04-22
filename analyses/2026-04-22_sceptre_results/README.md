# Create BEDPE files for SCEPTRE discovery results

python3 scripts/python/sceptre_to_bedpe.py \
  --input-file results/2026-04-22_sceptre_results/day0_grna20/sceptre_discovery_results.csv \
  --tss-file annotations/genes/gencode.v43.protein_coding.TSS500bp.bed \
  --output-file results/2026-04-22_sceptre_results/day0_grna20/sceptre_discovery_results_day0_grna20.bedpe \
  --tss-midpoint \
  --element-midpoint

python3 scripts/python/sceptre_to_bedpe.py \
  --input-file results/2026-04-22_sceptre_results/day2_grna20/sceptre_discovery_results.csv \
  --tss-file annotations/genes/gencode.v43.protein_coding.TSS500bp.bed \
  --output-file results/2026-04-22_sceptre_results/day2_grna20/sceptre_discovery_results_day2_grna20.bedpe \
  --tss-midpoint \
  --element-midpoint

python3 scripts/python/sceptre_to_bedpe.py \
  --input-file results/2026-04-22_sceptre_results/day4_grna20/sceptre_discovery_results.csv \
  --tss-file annotations/genes/gencode.v43.protein_coding.TSS500bp.bed \
  --output-file results/2026-04-22_sceptre_results/day4_grna20/sceptre_discovery_results_day4_grna20.bedpe \
  --tss-midpoint \
  --element-midpoint



## Annotate SCEPTRE discovery results
python3 scripts/python/sceptre_pairs_annotations.py \
    --input-file results/2026-04-22_sceptre_results/day0_grna20/sceptre_discovery_results.csv \
    --tss-file annotations/genes/gencode.v43.protein_coding.TSS500bp.bed \
    --element-annotations results/03-2025-10-24_create_sceptre_input/elements.annotated.provenance.with_tss_distance.tsv \
    --quadrant-annotations results/2026-04-21_load_counts_and_qc/day0_ribo20_tumi55_umi1500/gene_umi_comparison_quadrants.tsv \
    --gene-trends-file results/gene_classification_by_Noam/gene_grouping_final.tsv \
    --output-file results/2026-04-22_sceptre_results/annotated_results/sceptre_discovery_results_day0_grna20_annotated.tsv \
    --tss-midpoint \
    --element-midpoint

python3 scripts/python/sceptre_pairs_annotations.py \
    --input-file results/2026-04-22_sceptre_results/day2_grna20/sceptre_discovery_results.csv \
    --tss-file annotations/genes/gencode.v43.protein_coding.TSS500bp.bed \
    --element-annotations results/03-2025-10-24_create_sceptre_input/elements.annotated.provenance.with_tss_distance.tsv \
    --quadrant-annotations results/2026-04-21_load_counts_and_qc/day2_ribo16_tumi45_umi1000/gene_umi_comparison_quadrants.tsv \
    --gene-trends-file results/gene_classification_by_Noam/gene_grouping_final.tsv \
    --output-file results/2026-04-22_sceptre_results/annotated_results/sceptre_discovery_results_day2_grna20_annotated.tsv \
    --tss-midpoint \
    --element-midpoint

python3 scripts/python/sceptre_pairs_annotations.py \
    --input-file results/2026-04-22_sceptre_results/day4_grna20/sceptre_discovery_results.csv \
    --tss-file annotations/genes/gencode.v43.protein_coding.TSS500bp.bed \
    --element-annotations results/03-2025-10-24_create_sceptre_input/elements.annotated.provenance.with_tss_distance.tsv \
    --quadrant-annotations results/2026-04-21_load_counts_and_qc/day4_ribo6_tumi25_umi1000/gene_umi_comparison_quadrants.tsv \
    --gene-trends-file results/gene_classification_by_Noam/gene_grouping_final.tsv \
    --output-file results/2026-04-22_sceptre_results/annotated_results/sceptre_discovery_results_day4_grna20_annotated.tsv \
    --tss-midpoint \
    --element-midpoint