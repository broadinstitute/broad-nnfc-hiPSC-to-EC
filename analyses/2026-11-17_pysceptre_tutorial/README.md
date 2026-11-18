# pysceptre tutorial

A marimo notebook, `pysceptre-tutorial-engreitzlab.py`, that analyses day 0 of
the hiPSC-to-EC screen with pysceptre. It covers:

- RNA and guide QC;
- guide assignment;
- the positive-control, host-gene and discovery pairs;
- the calibration and discovery tests;
- power;
- whether long-range links clear the trans level.

## Install

Install [pixi](https://pixi.sh), then run this from this folder:

```sh
pixi install
```

This builds the environment from `pixi.lock`, which is locked for macOS on Apple
silicon and Linux x86-64. Python 3.13 comes from conda-forge; pysceptre 0.6.0
and the other packages come from PyPI, at the versions the notebook's own
header pins.

## Run

```sh
pixi run edit
```

This opens the notebook in the browser. Extra marimo flags go after the task,
for example `pixi run edit --headless --port 2720` on a remote machine.

The first run downloads `day0_pass_qc.h5mu` (1.0 GB) into `data/` and checks
its md5. Later runs read that local copy. `data/` is gitignored.

You can also run the notebook without pixi, using
`uvx marimo edit --sandbox pysceptre-tutorial-engreitzlab.py`. It reads the
same pins from the notebook's PEP 723 header, so change `pixi.toml` and that
header together.

## Settings

- `N_CELLS = 50_000`, at the top: a random sample of cells (`SAMPLE_SEED = 0`),
  so a training run is quick. With `None` the notebook analyses every cell,
  and only then do the numbers match the full screen.
- `N_JOBS = 8`, in the calibration cell: the CPUs pysceptre uses. Lower it on a
  smaller machine.
- `TRANS_CHUNKS_PER_CALL = 4`, in the long-range cell: how many of
  discovery's fit chunks each trans call fits. Fewer holds less in memory at
  once.

## Where the data comes from

`slim_day0.py` built `day0_pass_qc.h5mu` from
`results/2026-11-12_day_objects/days/day0.h5mu`. It keeps the cells that pass
QC and the total counts. The gene `var` gains `gene_tss` and `gene_strand` from
GENCODE v43. The script records how the file was made; it is not a step of the
tutorial. It needs crispr-analysis-utils and the local day object.
