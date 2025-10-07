#!/usr/bin/env bash
set -euo pipefail

OUTDIR="../../data/2025-09-02"
if [[ -d $OUTDIR ]]; then
  echo "Directory already exists"
else
  mkdir -p $OUTDIR
fi

#Get the list of paths to the h5 files on GCP.
gcloud storage ls "gs://250806_d0_250617_ec-large-scale/cellranger-outputs/D0-outputs/**/per_sample_outs/**/count/sample_filtered_feature_bc_matrix.h5" > h5-uris.txt

# Paste your URIs between EOFs
while IFS= read -r uri; do
  # Extract sample ID like 150-D0-r2-A6 from the path
  sample_id="$(grep -oE '[0-9]+-D0-r2-[AB][0-9]+' <<< "$uri" | head -n1)"
  dest="${OUTDIR}/${sample_id}_sample_filtered_feature_bc_matrix.h5"
  echo "Downloading $uri -> $dest"
  gcloud storage cp "$uri" "$dest"
done < h5-uris.txt


