# QC of the day objects

One marimo notebook per day, starting from that day's object in
`results/2026-11-12_day_objects/days/`. Each notebook:

- keeps the cells `flag_cells.py` marked `pass_qc` that also pass the R QC's
  three limits for that day (`analyses/2026-01-20_load_counts_and_qc`): gene
  UMIs, readout share and ribosomal share (`obs["keep"]`);
- adds `gene_tss` and `gene_strand` to the genes, from the GENCODE v43 store;
- runs the pysceptre tutorial's steps (`analyses/2026-11-17_pysceptre_tutorial`)
  through RNA QC, guide QC, guide assignment, the pairs and the covariates;
- writes the cells kept to `results/2026-11-23_day_objects_qc/day{N}_pass_qc.h5mu`
  when its "Write the h5mu" button is pressed.

The calibration check, power and discovery are not here: they come later, in
their own analysis folder.

| Notebook     | Day | Gene UMIs | Readout share | Ribosomal share |
|--------------|-----|-----------|---------------|-----------------|
| `day0_qc.py` | 0   | ≥ 1,500   | ≥ 55%         | ≥ 20%           |
| `day2_qc.py` | 2   | ≥ 1,000   | ≥ 45%         | ≥ 16%           |
| `day4_qc.py` | 4   | ≥ 1,000   | ≥ 25%         | ≥ 6%            |

The code is duplicated across days on purpose, so each day's thresholds can
differ.

## Install

Install [pixi](https://pixi.sh), then run this from this folder:

```sh
pixi install
```

Python 3.13 comes from conda-forge; pysceptre 0.6.0, crispr-analysis-utils
(from git, at a pinned commit) and the other packages come from PyPI, at the
versions the notebooks' headers pin.

## Run

```sh
pixi run edit day0_qc.py
```

The notebook reads the day object and
`~/Annotations/HG38/GENCODE/gencode.v43.annotation.parquet` from this machine;
nothing is downloaded.

You can also run it without pixi, with `uvx marimo edit --sandbox day0_qc.py`,
which reads the pins from the notebook's PEP 723 header. Change `pixi.toml` and
that header together.

## Settings

- `N_CELLS`, at the top: `None` analyses every cell kept; a number analyses a
  random sample of them (`SAMPLE_SEED = 0`), for a quick run.
- `MIN_GENE_UMIS`, `MIN_READOUT_PCT`, `MIN_RIBO_PCT`, in "Keep the cells that
  pass QC": the three limits, per day.
