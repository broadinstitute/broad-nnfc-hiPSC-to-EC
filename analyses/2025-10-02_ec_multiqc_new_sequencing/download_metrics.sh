#!/usr/bin/env bash
set -euo pipefail

# ---- EDIT THESE ----
BUCKET="gs://fc-76565551-28c3-4b6d-a048-e272103bcbd1/cellranger_scatter_outputs"
SAMPLES_FILE="samples.txt"
OUTDIR="../../results/2025-10-02_ec_multiqc_new_sequencing/metrics_qc"

while read -r sample || [[ -n "${sample:-}" ]]; do
  [[ -z "$sample" ]] && continue
  src="${BUCKET}/${sample}/outs/metrics_summary.csv"
  dst="${OUTDIR}/${sample}_metrics_summary.csv"

  echo ">> ${sample}"
  gcloud storage cp "$src" "$dst" \
  || echo "WARN: missing or inaccessible: $src"
done < "$SAMPLES_FILE"
