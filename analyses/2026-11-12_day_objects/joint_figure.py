"""Gene UMIs against guide UMIs per cell, one density figure per day.

Reads only the cell tables of ``days/day<N>.h5mu`` (not the matrices) and
publishes ``figures/joint_day<N>`` with `cau.qc.metric_density`. A cell the
guide modality does not hold has 0 guide UMIs.

Usage::

    pixi run --manifest-path ~/GitHub/crispr-analysis-utils/pyproject.toml -e dev \
        python joint_figure.py [--days 0 2 4] [--gene-threshold N] [--png DIR]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import h5py
import pandas as pd
from anndata.io import read_elem

from crispr_analysis_utils import qc

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results" / "2026-11-12_day_objects"


def cells(path: Path) -> pd.DataFrame:
    """Gene and guide UMIs per cell, from the day object's obs tables."""
    with h5py.File(path, "r") as f:
        gene = read_elem(f["mod/gene/obs"])
        guide = read_elem(f["mod/guide/obs"])
    return pd.DataFrame(
        {
            "gene_umis": gene["total_counts"],
            "guide_umis": guide["total_guide_umis"].reindex(gene.index).fillna(0),
        }
    )


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--days", nargs="*", type=int, default=[0, 2, 4])
    parser.add_argument("--gene-threshold", type=float)
    parser.add_argument("--png", type=Path)
    args = parser.parse_args(argv)

    for day in args.days:
        frame = cells(OUT / "days" / f"day{day}.h5mu")
        figure, numbers, caption = qc.metric_density(
            frame,
            "gene_umis",
            "guide_umis",
            x_threshold=args.gene_threshold,
            xlabel="gene UMIs per cell",
            ylabel="guide UMIs per cell",
        )
        caption = f"Day {day}. " + caption
        qc.publish(OUT / "figures" / f"joint_day{day}", figure, numbers, caption)
        if args.png:
            args.png.mkdir(parents=True, exist_ok=True)
            figure.savefig(args.png / f"joint_day{day}.png", dpi=130, bbox_inches="tight")
        print(caption)
    return 0


if __name__ == "__main__":
    sys.exit(main())
