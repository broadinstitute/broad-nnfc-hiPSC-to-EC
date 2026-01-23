#!/usr/bin/env bash
set -euo pipefail


# Error handler to catch any failure
error_handler() {
  local exit_code=$?
  echo "ERROR: Script failed at line $1 with exit code $exit_code" | tee -a "$LOG"

  # Remove corrupted bigWig if present
  if [[ -n "${OUT_BW:-}" && -f "$OUT_BW" ]]; then
    echo "Removing corrupted output file: $OUT_BW" | tee -a "$LOG"
    rm -f "$OUT_BW"
  fi
    
  # Remove temporary bedGraph if present
  if [[ -n "${TMP_BG:-}" && -f "$TMP_BG" ]]; then
    echo "Removing temporary bedGraph file: $TMP_BG" | tee -a "$LOG"
    rm -f "$TMP_BG"
  fi
  if [[ -n "${TMP_BG_FILT:-}" && -f "$TMP_BG_FILT" ]]; then
    echo "Removing filtered temporary bedGraph file: $TMP_BG_FILT" | tee -a "$LOG"
    rm -f "$TMP_BG_FILT"
  fi


  echo "Check the log file ($LOG) for details." | tee -a "$LOG"
  exit $exit_code
}
trap 'error_handler $LINENO' ERR


# Author: Eugenio Mattei
# Affiliation: Broad Institute of MIT and Harvard
# Date: 2025-11-09
# License: MIT
# Requirements:
#   - samtools
#   - bedtools
#   - bedGraphToBigWig (from UCSC Kent utils)
# Description:
#   Create a CPM-normalized bigWig file from a BAM file.
#   CPM normalization scales the coverage such that the total number of
#   mapped reads equals 1 million.
# Usage: bigwig_cpm.sh input.bam genome.sizes
# Output:
#   - <prefix>_CPM.bw
#   - <prefix>.log

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 input.bam genome.sizes" >&2
  exit 1
fi

BAM="$1"
GENOME_SIZES="$2"

if [[ ! -f "$BAM" ]]; then
  echo "Error: BAM file '$BAM' not found" >&2
  exit 1
fi

if [[ ! -f "$GENOME_SIZES" ]]; then
  echo "Error: genome sizes file '$GENOME_SIZES' not found" >&2
  exit 1
fi

# Prefix from BAM file name (strip path and .bam)
PREFIX="$(basename "${BAM%.bam}")"
LOG="${PREFIX}.log"
TMP_BG=$(mktemp --suffix=".bedgraph")
OUT_BW="${PREFIX}_CPM.bw"

# Archive existing log if present
if [[ -f "$LOG" ]]; then
  TS=$(date +%Y%m%d_%H%M%S)
  mv "$LOG" "${PREFIX}_${TS}.log"
fi

# Start fresh log file
: > "$LOG"


# Start log
{
  echo "[$(date)] Starting CPM bigWig generation"
  echo "BAM:          $BAM"
  echo "Genome sizes: $GENOME_SIZES"
  echo "Prefix:       $PREFIX"
  echo
} >> "$LOG"

echo "Computing total mapped reads..." | tee -a "$LOG"

# Total mapped reads from samtools idxstats (3rd column)
TOTAL_READS=$(samtools idxstats "$BAM" 2>>"$LOG" | \
              awk '{sum += $3} END {print sum}')

if [[ -z "$TOTAL_READS" || "$TOTAL_READS" -eq 0 ]]; then
  echo "Error: total mapped reads is zero or undefined. Check BAM/index." | tee -a "$LOG" >&2
  exit 1
fi

echo "Computing scaling factor..." | tee -a "$LOG"

# Compute CPM scale factor: 1e6 / total_mapped_reads
SCALE=$(awk -v t="$TOTAL_READS" 'BEGIN {printf "%.10f", 1000000 / t}')

{
  echo "Total mapped reads: $TOTAL_READS"
  echo "CPM scale factor:   $SCALE"
  echo
  echo "Running: bedtools genomecov -> temp bedGraph -> filter -> bedGraphToBigWig"
} >> "$LOG"

# temp files
TMP_BG=$(mktemp --suffix=".bedgraph")
TMP_BG_FILT=$(mktemp --suffix=".bedgraph")

{
  echo "Temporary bedGraph (raw):   $TMP_BG"
  echo "Temporary bedGraph (filt):  $TMP_BG_FILT"
} >> "$LOG"

# 1) CPM coverage
echo "Computing CPM coverage with bedtools genomecov..." | tee -a "$LOG"

bedtools genomecov -ibam "$BAM" -bg -scale "$SCALE" 2>>"$LOG" > "$TMP_BG"

# 2) filter bedGraph to only chromosomes present in genome.sizes
echo "Filtering bedGraph to chromosomes in genome sizes..." | tee -a "$LOG"
awk '
  NR==FNR {allowed[$1]=1; next}
  ($1 in allowed)
' "$GENOME_SIZES" "$TMP_BG" > "$TMP_BG_FILT"

# 3) convert to bigWig
echo "Converting filtered bedGraph to bigWig..." | tee -a "$LOG"

bedGraphToBigWig "$TMP_BG_FILT" "$GENOME_SIZES" "$OUT_BW" 2>>"$LOG"

# cleanup on success
rm -f "$TMP_BG" "$TMP_BG_FILT"

{
  echo
  echo "[$(date)] Done."
  echo "Output bigWig: $OUT_BW"
} >> "$LOG"

echo "Finished. bigWig: $OUT_BW  |  Log: $LOG"
