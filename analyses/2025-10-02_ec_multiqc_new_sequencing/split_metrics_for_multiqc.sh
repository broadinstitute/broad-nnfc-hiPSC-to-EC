# split_metrics_for_multiqc.sh
set -euo pipefail
IN="../../results/2025-10-02_ec_multiqc_new_sequencing/metrics_qc/metrics_qc_custom_mqc.tsv"
OUTDIR="../../results/2025-10-02_ec_multiqc_new_sequencing/mqc_plots"
mkdir -p "$OUTDIR"

python3 - "$IN" "$OUTDIR" <<'PY'
import csv, os, re, sys
src, outdir = sys.argv[1], sys.argv[2]

def slug(s):
    s = re.sub(r'[:()/]', '', s)
    s = re.sub(r'\s+', '_', s.strip())
    s = re.sub(r'[^A-Za-z0-9_.-]+', '_', s)
    return s[:80]

# read TSV (skip MultiQC header lines starting with '#')
with open(src, newline='') as f:
    lines = [l for l in f if not l.startswith('#')]
header = next(csv.reader([lines[0]], delimiter='\t'))
cols = header[1:]                # first col is 'Sample'
rows = list(csv.reader(lines[1:], delimiter='\t'))

samples = [r[0] for r in rows]

for idx, colname in enumerate(cols, start=1):
    slugname = slug(colname)
    fn = os.path.join(outdir, f"{slugname}_mqc.tsv")
    with open(fn, 'w', newline='') as fo:
        # MultiQC file header -> one plot per file
        fo.write(f'# id: "cellranger_metric_{slugname}"\n')
        fo.write(f'# section_name: "Cell Ranger: {colname}"\n')
        fo.write('# format: "tsv"\n')
        fo.write('# plot_type: "bargraph"\n')
        fo.write('# pconfig:\n')
        fo.write(f'#   ylab: "{colname}"\n')
        fo.write('#   col1_header: "Sample"\n')
        # data
        fo.write('Sample\tvalue\n')
        for r in rows:
            val = r[idx] if idx < len(r) else ''
            fo.write(f'{r[0]}\t{val}\n')
print(f"Wrote {len(cols)} plot files to {outdir}")
PY

