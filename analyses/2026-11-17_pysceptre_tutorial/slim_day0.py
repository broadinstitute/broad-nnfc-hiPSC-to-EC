"""Day 0 for the pysceptre tutorial: the cells that pass QC, total counts only.

Reads ``results/2026-11-12_day_objects/days/day0.h5mu`` and writes
``results/2026-11-17_pysceptre_tutorial/day0_pass_qc.h5mu``:

- cells: those with ``obs['pass_qc']``;
- ``gene``: X (total UMIs), obs and var, without the three kb layers. var
  gains ``gene_tss`` (0-based) and ``gene_strand`` from the GENCODE v43 store
  the v2 pairs were built with, so pairs can be rebuilt without a GTF;
- ``guide``: X (UMIs), the ``guide_assignment`` layer, obs, var and uns;
- ``uns``: the day object's, plus ``source``.

Usage::

    pixi run --manifest-path ~/GitHub/crispr-analysis-utils/pyproject.toml -e dev \
        python slim_day0.py
"""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import anndata as ad
import h5py
import mudata as md
import numpy as np
import pandas as pd
from anndata.io import read_elem, sparse_dataset
from crispr_analysis_utils.guide_annotation.db import load_store
from crispr_analysis_utils.pairs import gene_tss

REPO = Path(__file__).resolve().parents[2]
SOURCE = REPO / "results" / "2026-11-12_day_objects" / "days" / "day0.h5mu"
STORE = Path.home() / "Annotations" / "HG38" / "GENCODE" / "gencode.v43.annotation.parquet"
OUT = REPO / "results" / "2026-11-17_pysceptre_tutorial" / "day0_pass_qc.h5mu"


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(16 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def modality(f: h5py.File, name: str, cells: pd.Index, layers: tuple[str, ...] = ()) -> ad.AnnData:
    """One modality, restricted to `cells` and to the layers named."""
    group = f["mod"][name]
    obs = read_elem(group["obs"])
    rows = np.flatnonzero(obs.index.isin(cells))
    if rows.size != len(cells):
        raise SystemExit(f"{name}: {len(cells) - rows.size:,} kept cells are not in this modality")
    return ad.AnnData(
        X=sparse_dataset(group["X"])[rows],
        obs=obs.iloc[rows],
        var=read_elem(group["var"]),
        layers={layer: sparse_dataset(group["layers"][layer])[rows] for layer in layers},
        uns=read_elem(group["uns"]),
    )


def main() -> int:
    with h5py.File(SOURCE, "r") as f:
        obs = read_elem(f["obs"])
        cells = obs.index[obs["pass_qc"].to_numpy(dtype=bool)]
        gene = modality(f, "gene", cells)
        guide = modality(f, "guide", cells, layers=("guide_assignment",))
        uns = read_elem(f["uns"])

    tss = gene_tss(load_store(STORE))
    ids = gene.var_names.str.replace(r"\.\d+$", "", regex=True)
    gene.var["gene_tss"] = pd.array(tss["tss"].reindex(ids), dtype="Int64")
    gene.var["gene_strand"] = pd.Categorical(tss["strand"].reindex(ids))
    no_tss = int(gene.var["gene_tss"].isna().sum())
    no_tss_readout = int((gene.var["gene_tss"].isna() & gene.var["readout"]).sum())
    if no_tss_readout:
        raise SystemExit(f"{no_tss_readout} readout genes have no TSS in {STORE.name}")

    data = md.MuData({"gene": gene, "guide": guide})
    data.obs = obs.reindex(data.obs_names)
    data.uns.update(uns)
    data.uns["source"] = {
        "file": str(SOURCE.relative_to(REPO)),
        "md5": md5(SOURCE),
        "cells": "obs['pass_qc']",
        "store": STORE.name,
        "store_md5": md5(STORE),
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    partial = OUT.with_suffix(".partial.h5mu")
    data.write_h5mu(partial, compression="gzip")
    partial.rename(OUT)
    print(
        f"{data.n_obs:,} of {len(obs):,} cells, {gene.n_vars:,} genes ({no_tss:,} without a TSS,"
        f" none in the readout), {guide.n_vars:,} guides, {gene.X.nnz:,} gene and"
        f" {guide.X.nnz:,} guide nonzeros, {OUT.stat().st_size / 1e9:.2f} GB -> {OUT}"
    )
    print(f"md5 {md5(OUT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
