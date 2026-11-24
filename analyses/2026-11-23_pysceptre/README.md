# pysceptre on the day objects

One marimo notebook per day, starting from the cells that day's QC notebook
kept: `results/2026-11-23_day_objects_qc/day{N}_pass_qc.h5mu`, written by
`analyses/2026-11-23_day_objects_qc`. Each notebook:

- rebuilds the guide calls, the elements, the pairs, the pairwise QC and the
  covariates with the QC notebook's own cells, so they match it;
- runs the pysceptre tutorial's steps (`analyses/2026-11-17_pysceptre_tutorial`):
  the calibration check, the power check, discovery, power and the long-range
  check;
- corrects discovery with Benjamini-Hochberg and with IHW (pyIHW, weighted by
  distance to the TSS, at 20 and 5 strata), compares the three at the trans
  level, and runs the DC-TAP paper's direct/indirect model.

Six cells take long and start disabled: `calibration_nt`, `calibration_safe`,
`power_check`, `discovery_run`, `gene_fits` and `trans_results`. Enable them in
that order. Nothing is written to `results/` yet.

Days 2 and 4 are day 0's notebook with that day's file and labels. They keep
day 0's pairwise QC cuts and IHW settings, so the days compare.

| Notebook            | Day |
|---------------------|-----|
| `day0_pysceptre.py` | 0   |
| `day2_pysceptre.py` | 2   |
| `day4_pysceptre.py` | 4   |

## Install and run

Install [pixi](https://pixi.sh), then run this from this folder:

```sh
pixi install
pixi run edit day0_pysceptre.py
```

Or, without pixi, `uvx marimo edit --sandbox day0_pysceptre.py`, which reads the
pins from the notebook's PEP 723 header. Change `pixi.toml` and that header
together.
