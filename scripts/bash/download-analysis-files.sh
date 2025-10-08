#!/usr/bin/env bash
set -euo pipefail

#OUTDIR="${OUTDIR:-out}"
OUTDIR="../../data/h5"
mkdir -p "$OUTDIR"

while IFS= read -r uri; do
  [[ -z "$uri" ]] && continue

  # Strip the trailing file and pull the parent directory name
  sample_dir="${uri%/outs/filtered_feature_bc_matrix.h5}"
  sample_id="${sample_dir##*/}"

  if [[ -z "${sample_id}" || "$sample_dir" == "$uri" ]]; then
    echo "⚠️  Could not parse sample_id from: $uri" >&2
    continue
  fi

  dest="${OUTDIR}/${sample_id}_sample_filtered_feature_bc_matrix.h5"
  echo "Downloading $uri -> $dest"
  gcloud storage cp "$uri" "$dest"
done < h5_uris.txt
