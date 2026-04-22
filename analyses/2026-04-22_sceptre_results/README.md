

python3 scripts/python/sceptre_to_bedpe.py \
  --input-file results/2026-04-22_sceptre_results/day0_grna20/sceptre_discovery_results.csv \
  --tss-file annotations/genes/gencode.v43.protein_coding.TSS500bp.bed \
  --genome-file annotations/GRCh38_EBV.chrom.sizes.no.alt.tsv \
  --output-file results/2026-04-22_sceptre_results/day0_grna20/sceptre_discovery_results_day0_grna20.bedpe \
  --tss-midpoint \
  --element-midpoint

python3 scripts/python/sceptre_to_bedpe.py \
  --input-file results/2026-04-22_sceptre_results/day2_grna20/sceptre_discovery_results.csv \
  --tss-file annotations/genes/gencode.v43.protein_coding.TSS500bp.bed \
  --genome-file annotations/GRCh38_EBV.chrom.sizes.no.alt.tsv \
  --output-file results/2026-04-22_sceptre_results/day2_grna20/sceptre_discovery_results_day2_grna20.bedpe \
  --tss-midpoint \
  --element-midpoint

python3 scripts/python/sceptre_to_bedpe.py \
  --input-file results/2026-04-22_sceptre_results/day4_grna20/sceptre_discovery_results.csv \
  --tss-file annotations/genes/gencode.v43.protein_coding.TSS500bp.bed \
  --genome-file annotations/GRCh38_EBV.chrom.sizes.no.alt.tsv \
  --output-file results/2026-04-22_sceptre_results/day4_grna20/sceptre_discovery_results_day4_grna20.bedpe \
  --tss-midpoint \
  --element-midpoint