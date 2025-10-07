#!/usr/bin/env bash
set -euo pipefail

INDIR="../../results/2025-10-02_ec_multiqc_new_sequencing/metrics_qc"
OUT="../../results/2025-10-02_ec_multiqc_new_sequencing/metrics_qc/metrics_qc_custom_mqc.tsv"

# ---- MultiQC header ----
cat >"$OUT" <<'EOF'
# id: "cellranger_metrics_summary"
# section_name: "Cell Ranger metrics_summary.csv"
# description: "Parsed metrics_summary per sample (cleaned numbers: no thousands commas, no %)"
# format: "tsv"
# plot_type: "table"
# pconfig:
#   col1_header: "Sample"
EOF

# POSIX glob expansion; error if no matches
set -- "$INDIR"/*_metrics_summary.csv
[ "$1" = "$INDIR/*_metrics_summary.csv" ] && { echo "No *_metrics_summary.csv in $INDIR" >&2; exit 1; }

# Feed the list (sorted) to Python for proper CSV parsing
printf '%s\n' "$@" | sort | python3 - "$OUT" <<'PY'
import csv, sys, os, re

out_path = sys.argv[1]
files = [line.strip() for line in sys.stdin if line.strip()]
if not files:
    sys.exit("No input files from stdin")

def clean_cell(v: str) -> str:
    v = v.strip().replace("\r","")
    if v.endswith("%"):
        v = v[:-1]
    # remove thousands separators only for numeric tokens like 1,234 or 1,234.56
    if re.fullmatch(r"\d[\d,]*(\.\d+)?", v):
        v = v.replace(",", "")
    return v

# Read header from first file (CSV reader handles commas & quoting)
with open(files[0], newline="") as fh:
    header = next(csv.reader(fh))

# Write header
with open(out_path, "a", newline="") as oh:
    oh.write("Sample\t" + "\t".join(header) + "\n")

# Process rows
expected = len(header)
with open(out_path, "a", newline="") as oh:
    for f in files:
        with open(f, newline="") as fh:
            rdr = csv.reader(fh)
            hdr = next(rdr, None)
            row = next(rdr, None)
        if hdr is None or row is None:
            sys.exit(f"{f}: missing header or data row")
        if len(hdr) != expected:
            sys.exit(f"Header mismatch in {f}: got {len(hdr)} cols, expected {expected}")
        if len(row) != expected:
            # Show the raw row joined by commas to help debugging
            sys.exit(f"Row colcount mismatch in {f}: got {len(row)}, expected {expected}\nRow was: {','.join(row)}")

        row = [clean_cell(x) for x in row]
        sample = os.path.basename(f).replace("_metrics_summary.csv","")
        oh.write(sample + "\t" + "\t".join(row) + "\n")

print(f"Wrote {out_path}")
PY

# Quick integrity check: all lines have same #fields
awk -F'\t' '!/^#/ {print NF}' "$OUT" | sort -u
