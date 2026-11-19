"""Flag every cell of each day object for RNA and guide quality; drop nothing.

Writes, in place, into ``days/day<N>.h5mu``:

- root ``obs``:
    - ``rna_umi_cut``: gene UMIs at or above the day's cut;
    - ``pass_rna``: past the cut and inside every MAD band;
    - ``pass_guide``: at least one guide assigned (the 20-UMI pysceptre calls);
    - ``pass_qc``: both. Inference reads the cells where it is true.
- ``mod['gene'].obs``: one nullable-boolean column per rule, and
  ``pct_counts_ribo_rmd``. ``NA`` where no band was measured: below the cut,
  or in a channel whose MAD for that metric is 0, which no cell fails.
- ``uns['qc_rules']``: every rule of every channel, with its bounds and
  the cells it failed; ``uns['qc_flags']``: the settings.

The MAD bands are measured per channel, on the cells past the cut: a band
measured with the ambient barcodes in it describes the ambient barcodes.

The R QC's per-day limits on the targeted and ribosomal fractions
(``analyses/2026-01-20_load_counts_and_qc``) are measured too, as
``*_limit`` columns, but are not part of ``pass_rna``.

Also writes ``cascade.tsv``: per channel, the cells left after each step.

Usage::

    pixi run --manifest-path ~/GitHub/crispr-analysis-utils/pyproject.toml -e dev \
        python flag_cells.py [--days 0 2 4]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata.io import read_elem, write_elem

from crispr_analysis_utils import qc

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results" / "2026-11-12_day_objects"

GENE_UMI_CUT = {0: 1000, 2: 1000, 4: 700}
NMADS = 5
BAND = ("total_counts", "n_genes_by_counts", "pct_counts_mito")
# The R QC's thresholds, read off each day's distributions; ribosomal on its
# own wider definition. Measured for comparison, not applied.
R_LIMITS = {
    0: {"pct_counts_targeted": 55, "pct_counts_ribo_rmd": 20},
    2: {"pct_counts_targeted": 45, "pct_counts_ribo_rmd": 16},
    4: {"pct_counts_targeted": 25, "pct_counts_ribo_rmd": 6},
}
RIBO_RMD = r"^RPS|^RPL|^RRP|^MRP|^FAU"
METRICS = BAND + ("pct_counts_targeted", "pct_counts_ribo_rmd")


def ribo_rmd(f: h5py.File, gene_obs: pd.DataFrame) -> np.ndarray:
    """Percent of each cell's gene UMIs on the R QC's ribosomal genes."""
    var = read_elem(f["mod/gene/var"])
    mask = var["symbol"].astype("string").fillna("").str.contains(RIBO_RMD, regex=True)
    x = f["mod/gene/X"]
    matrix = sp.csr_matrix(
        (x["data"][:], x["indices"][:], x["indptr"][:]), shape=tuple(x.attrs["shape"])
    )
    ribo = np.asarray(matrix @ mask.to_numpy().astype(np.float64)).ravel()
    return 100 * ribo / gene_obs["total_counts"].to_numpy()


def replace(group: h5py.Group, key: str, value) -> None:
    """Write `value` at `key`, swapping it in only once it is written."""
    write_elem(group, f"{key}__new", value)
    if key in group:
        del group[key]
    group.move(f"{key}__new", key)


def flag(day: int) -> pd.DataFrame:
    path = OUT / "days" / f"day{day}.h5mu"
    cut = GENE_UMI_CUT[day]
    with h5py.File(path, "r+") as f:
        root = read_elem(f["obs"])
        gene = read_elem(f["mod/gene/obs"])
        guide = read_elem(f["mod/guide/obs"])
        if not gene.index.equals(root.index):
            raise SystemExit(f"day {day}: gene and root obs differ")
        gene["pct_counts_ribo_rmd"] = ribo_rmd(f, gene)

        above = gene["total_counts"] >= cut
        rules_all, columns = [], {}
        for ms, cells in gene[above].groupby("measurement_set", observed=True):
            passes, rules = qc.cell_filters(
                cells[list(METRICS)],
                nmads=NMADS,
                band=BAND,
                limits={m: (low, None) for m, low in R_LIMITS[day].items()},
            )
            rules_all.append(rules.assign(measurement_set=ms, n_cells=len(cells)))
            for name in passes:
                columns.setdefault(name, []).append(passes[name])
        for name, parts in columns.items():
            gene[name] = pd.concat(parts).reindex(gene.index).astype("boolean")

        mad = [f"{m}_mad" for m in BAND]
        limit = [f"{m}_limit" for m in R_LIMITS[day]]
        assigned = guide["n_guides_assigned"].reindex(root.index).fillna(0)
        root["rna_umi_cut"] = above.to_numpy()
        # NA past the cut: the channel has no band for that metric (a MAD of 0).
        root["pass_rna"] = (above & gene[mad].fillna(True).all(axis=1)).to_numpy()
        root["pass_guide"] = (assigned >= 1).to_numpy()
        root["pass_qc"] = root["pass_rna"] & root["pass_guide"]

        replace(f, "obs", root)
        replace(f["mod/gene"], "obs", gene)
        replace(f["uns"], "qc_rules", pd.concat(rules_all, ignore_index=True))
        replace(
            f["uns"],
            "qc_flags",
            {
                "gene_umi_cut": cut,
                "nmads": NMADS,
                "band": list(BAND),
                "band_measured_on": "each channel's cells at or above gene_umi_cut",
                "pass_rna": "rna_umi_cut and every *_mad rule",
                "pass_guide": "n_guides_assigned >= 1 (pysceptre thresholding, 20 UMIs)",
                "pass_qc": "pass_rna and pass_guide",
                "not_applied": "the *_limit rules: the R QC's per-day limits, for comparison",
            },
        )

    r_limits = gene[limit].fillna(False).all(axis=1)
    table = pd.DataFrame(
        {
            "barcodes": 1,
            "gene_umi_cut": above,
            "pass_rna": root["pass_rna"],
            "pass_qc": root["pass_qc"],
            "pass_qc_and_r_limits": root["pass_qc"] & r_limits,
            "pass_rna_no_guide": root["pass_rna"] & ~root["pass_guide"],
        },
        index=root.index,
    ).astype(int)
    for m in BAND:
        table[f"fail_{m}_mad"] = (above & ~gene[f"{m}_mad"].fillna(True)).astype(int)
    cascade = table.groupby(root["measurement_set"].to_numpy()).sum()
    cascade.index.name = "measurement_set"
    cascade = cascade.reset_index()
    cascade.insert(1, "sample", cascade["measurement_set"].map(
        root.groupby("measurement_set", observed=True)["sample"].first()))
    cascade.insert(1, "day", day)
    return cascade


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--days", nargs="*", type=int, default=sorted(GENE_UMI_CUT))
    args = parser.parse_args(argv)

    tables = []
    for day in args.days:
        cascade = flag(day)
        tables.append(cascade)
        total = cascade.drop(columns=["measurement_set", "sample", "day"]).sum()
        print(f"day {day} (cut {GENE_UMI_CUT[day]}):", total.to_dict())
    path = OUT / "cascade.tsv"
    old = pd.read_csv(path, sep="\t") if path.exists() else None
    new = pd.concat(tables, ignore_index=True)
    if old is not None:
        new = pd.concat([old[~old["day"].isin(args.days)], new], ignore_index=True)
    new.sort_values(["day", "sample"]).to_csv(path, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
