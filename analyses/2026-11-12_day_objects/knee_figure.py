"""The knees of every prepared channel, gene and guide, coloured by day.

Reads ``knee_curves/*.tsv.gz`` and ``channels.tsv`` from `prepare_channels.py`
and publishes ``figures/knee_<modality>`` (PDF, numbers, caption) with
`cau.qc.knee_overlay`. Run it again as channels finish.

Usage::

    pixi run --manifest-path ~/GitHub/crispr-analysis-utils/pyproject.toml -e dev \
        python knee_figure.py [--png DIR]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from crispr_analysis_utils import qc

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results" / "2026-11-12_day_objects"
FLOOR = 100


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--png", type=Path, help="Also write PNGs here, to look at.")
    args = parser.parse_args(argv)

    summary = pd.read_csv(OUT / "channels.tsv", sep="\t")
    day = dict(zip(summary["measurement_set"], summary["day"], strict=True))
    curves = pd.concat(
        [pd.read_csv(p, sep="\t") for p in sorted((OUT / "knee_curves").glob("*.tsv.gz"))],
        ignore_index=True,
    )
    curves = curves[curves["measurement_set"].isin(day)]
    curves["day"] = pd.Categorical(
        "day " + curves["measurement_set"].map(day).astype(str),
        categories=["day 0", "day 2", "day 4"],
    )
    curves = curves.rename(columns={"measurement_set": "channel"})

    for modality, label, floor in (
        ("gene", "gene UMIs per barcode", FLOOR),
        ("guide", "guide UMIs per barcode", None),
    ):
        figure, numbers, caption = qc.knee_overlay(
            curves[curves["modality"] == modality], floor=floor, ylabel=label
        )
        qc.publish(OUT / "figures" / f"knee_{modality}", figure, numbers, caption)
        if args.png:
            args.png.mkdir(parents=True, exist_ok=True)
            figure.savefig(args.png / f"knee_{modality}.png", dpi=150)
        print(f"knee_{modality}: {caption}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
