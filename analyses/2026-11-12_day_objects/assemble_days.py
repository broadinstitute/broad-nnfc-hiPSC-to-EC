"""One h5mu per day, from the channels `prepare_channels.py` wrote.

Each day's channels are joined cell-wise. The object keeps every gene and
guide, X and the three kb layers, the 20-UMI guide calls, and each cell's
channel metadata and QC metrics. Nothing is filtered beyond the 100-UMI
floor each channel was cut at; the QC rules are added later as flags.

``uns`` carries what a reader of the day needs without the channel files:

- ``channels``: one row per channel (`channels.tsv`);
- ``knee_curves``: every channel's barcode-rank curve, gene and guide, over
  all its barcodes, so the knee can be redrawn after the cut;
- ``qc``: each channel's floor record, by measurement set;
- ``reference``: the GTF and genome.

``mod['guide'].uns`` carries the guide annotation and element mapping
tables once, from the design folder.

Usage::

    pixi run --manifest-path ~/GitHub/crispr-analysis-utils/pyproject.toml -e dev \
        python assemble_days.py [--days 0 2 4] [--partial] [--out DIR]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results" / "2026-11-12_day_objects"
DESIGN = REPO / "results" / "2026-11-18_guide_design_v3"
EXPECTED = {0: 36, 2: 35, 4: 39}
CHANNEL_COLUMNS = ("replicate", "day", "chip", "batch", "channel", "sample", "channel_name")


def assemble(day: int, summary: pd.DataFrame):
    """Join one day's channel files into a MuData."""
    import anndata as ad
    import mudata as md

    rows = summary[summary["day"] == day].sort_values("measurement_set")
    genes, guides, obs, floors = [], [], [], {}
    reference = None
    for ms in rows["measurement_set"]:
        data = md.read_h5mu(OUT / "channels" / f"{ms}.h5mu")
        if genes and not data["gene"].var_names.equals(genes[0].var_names):
            raise SystemExit(f"{ms}: genes differ from the day's first channel")
        if guides and not data["guide"].var_names.equals(guides[0].var_names):
            raise SystemExit(f"{ms}: guides differ from the day's first channel")
        genes.append(data["gene"])
        guides.append(data["guide"])
        obs.append(
            data.obs[list(CHANNEL_COLUMNS)].assign(measurement_set=ms)
        )
        floors[ms] = dict(data.uns.get("qc", {}))
        reference = reference or dict(data.uns.get("reference", {}))
        assignment = dict(data["guide"].uns.get("guide_assignment", {}))

    gene = ad.concat(genes, merge="same")
    guide = ad.concat(guides, merge="same")
    guide.uns["guide_assignment"] = assignment
    guide.uns["annotation"] = pd.read_csv(
        DESIGN / "annotations" / "guide_gene_annotation.tsv", sep="\t", dtype=str
    )
    guide.uns["element_mapping"] = pd.read_csv(
        DESIGN / "elements" / "guide_element_mapping.tsv", sep="\t", dtype=str
    )

    data = md.MuData({"gene": gene, "guide": guide})
    root = pd.concat(obs)
    data.obs = root.reindex(data.obs_names)
    for column in data.obs:
        if data.obs[column].dtype == object:
            data.obs[column] = data.obs[column].astype("category")

    curves = pd.concat(
        [pd.read_csv(OUT / "knee_curves" / f"{ms}.tsv.gz", sep="\t") for ms in rows["measurement_set"]],
        ignore_index=True,
    )
    data.uns["channels"] = rows.reset_index(drop=True)
    data.uns["knee_curves"] = curves
    data.uns["qc"] = floors
    data.uns["reference"] = reference
    data.uns["day"] = day
    return data


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--days", nargs="*", type=int, default=sorted(EXPECTED))
    parser.add_argument("--partial", action="store_true", help="Allow a day with channels missing.")
    parser.add_argument("--out", type=Path, default=OUT / "days")
    args = parser.parse_args(argv)

    summary = pd.read_csv(OUT / "channels.tsv", sep="\t")
    for day in args.days:
        have = int((summary["day"] == day).sum())
        if have != EXPECTED[day] and not args.partial:
            raise SystemExit(f"day {day}: {have} of {EXPECTED[day]} channels prepared")
        data = assemble(day, summary)
        args.out.mkdir(parents=True, exist_ok=True)
        path = args.out / f"day{day}.h5mu"
        partial = path.with_suffix(".partial.h5mu")
        data.write_h5mu(partial, compression="gzip")
        partial.rename(path)
        held = np.asarray(data.obsm["guide"]).ravel()
        print(
            f"day {day}: {have} channels, {data.n_obs:,} cells, {int(held.sum()):,} with guide"
            f" counts, {path.stat().st_size / 1e9:.2f} GB -> {path}"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
