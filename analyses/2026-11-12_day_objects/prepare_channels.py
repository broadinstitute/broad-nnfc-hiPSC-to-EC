"""Every channel, built and cut to its candidate cells, one at a time.

For each of the 110 channels of the IGVF CRISPR pipeline run:

1. copy its kb RNA and guide counts from the staged run;
2. build the MuData with `cau mudata-builder`, on the v2 guide tables;
3. record the barcode-rank (knee) curve of both modalities over every
   barcode, so the knee figure sees the whole channel;
4. keep the barcodes with at least ``FLOOR`` gene UMIs, which drops empty
   droplets and no candidate cell (the floor sits below every knee; the
   knee figure is the check);
5. measure each kept cell (`qc.cell_fractions`, guide UMIs, guides seen) and
   call guides with pysceptre at ``THRESHOLD`` UMIs;
6. write ``channels/<measurement_set>.h5mu`` and delete the downloads.

Nothing is filtered beyond the floor: the QC rules become flags once the day
objects exist. X and the three kb layers are all kept. The guide annotation
and element mapping tables are the same for every channel (32 MB each time),
so the day objects carry them once instead.

Resumable: a channel whose file exists is skipped.

Usage::

    pixi run --manifest-path ~/GitHub/crispr-analysis-utils/pyproject.toml -e dev \
        python prepare_channels.py --workdir DIR [--channels IGVFDS... ...]
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "analyses" / "2026-10-15_counts_vs_cellranger"))
sys.path.insert(0, str(Path.home() / "GitHub" / "pysceptre" / "src"))

from compare_all_channels import SAMPLES, channels  # noqa: E402

STAGED = "gs://hgrm-seqera/igvf-crispr-all/per-channel"
DESIGN = REPO / "results" / "2026-11-18_guide_design_v3"
STORE = Path.home() / "Annotations/HG38/GENCODE/gencode.v43.annotation.parquet"
PANEL = REPO / "metadata" / "metadata_genes_final.csv"
GTF_ACCESSION = "IGVFFI9561BASO"
OUT = REPO / "results" / "2026-11-12_day_objects"
FLOOR = 100
THRESHOLD = 20
KNEE_POINTS = 400


def table() -> pd.DataFrame:
    """One row per channel, with the metadata columns of the sample sheet."""
    t = channels()
    meta = pd.read_csv(SAMPLES, sep="\t")[
        ["entity:endothelial_id", "replicate", "day", "chip", "batch", "channel"]
    ].drop_duplicates("entity:endothelial_id")
    t = t.drop(columns="day").merge(meta, on="entity:endothelial_id", how="left")
    if t[["replicate", "day", "chip", "batch", "channel"]].isna().any().any():
        raise SystemExit("a channel has no replicate, day, chip, batch or channel")
    return t.rename(columns={"entity:endothelial_id": "entity_id"})


def fetch(ms: str, workdir: Path) -> tuple[Path, Path]:
    """The channel's kb counts, copied from the staged run."""
    paths = []
    for kind, src in (
        ("rna", f"{STAGED}/{ms}/kb_rna/counts_unfiltered/adata.h5ad"),
        ("guide", f"{STAGED}/{ms}/kb_guide/counts_unfiltered_modified/adata.h5ad"),
    ):
        dst = workdir / ms / kind / "adata.h5ad"
        dst.parent.mkdir(parents=True, exist_ok=True)
        if not dst.exists():
            subprocess.run(["gcloud", "storage", "cp", "-q", src, str(dst)], check=True)
        paths.append(dst)
    return paths[0], paths[1]


def knee_curve(umis: np.ndarray, n_points: int = KNEE_POINTS) -> pd.DataFrame:
    """UMIs at log-spaced barcode ranks, over the barcodes with any UMI.

    Exact at every rank it keeps: the curve is sampled, not smoothed.
    """
    ranked = np.sort(umis[umis > 0])[::-1]
    if not len(ranked):
        return pd.DataFrame({"rank": [], "umis": []})
    ranks = np.unique(np.geomspace(1, len(ranked), n_points).round().astype(int))
    return pd.DataFrame({"rank": ranks, "umis": ranked[ranks - 1].astype(int)})


def narrow(x) -> sp.csr_matrix:
    """Counts as CSR in the narrowest unsigned type that holds them."""
    x = sp.csr_matrix(x)
    top = x.data.max() if x.nnz else 0
    return x.astype(np.uint16 if top <= np.iinfo(np.uint16).max else np.uint32)


def prepare(info, workdir: Path) -> dict:
    """Build, measure and cut one channel; return its summary row."""
    import mudata as md
    from pysceptre import assign_grnas_thresholding

    from crispr_analysis_utils import qc
    from crispr_analysis_utils.mudata_builder import run as build

    ms = info.measurement_sets
    rna, guide = fetch(ms, workdir)
    built = build(
        rna,
        guide,
        DESIGN / "igvf" / "guide_metadata.tsv",
        STORE,
        workdir / ms / "built",
        measurement_set=ms,
        gtf_accession=GTF_ACCESSION,
        readout_genes=PANEL,
        metadata={
            "replicate": str(int(info.replicate)),
            "day": str(int(info.day)),
            "chip": str(info.chip),
            "batch": str(int(info.batch)),
            "channel": str(int(info.channel)),
            "sample": str(info.sample),
            "channel_name": str(info.entity_id),
        },
    )
    data = built.data

    curves = []
    for modality in ("gene", "guide"):
        umis = qc.umis_per_cell(data, modality)
        curve = knee_curve(umis).assign(
            measurement_set=ms, modality=modality, n_barcodes=int((umis > 0).sum())
        )
        curves.append(curve)

    cut, _ = qc.filter_cells(data, FLOOR, modality="gene")
    gene = cut["gene"]
    fractions = qc.cell_fractions(cut, "gene").loc[gene.obs_names]
    for column in fractions:
        gene.obs[column] = fractions[column].to_numpy()

    guide_mod = cut["guide"]
    counts = sp.csr_matrix(guide_mod.X)
    called = assign_grnas_thresholding(
        counts.T.tocsr(), list(guide_mod.var_names), threshold=THRESHOLD
    )
    calls = sp.csr_matrix(called.assigned.T).astype(np.uint8)
    guide_mod.layers["guide_assignment"] = calls
    guide_mod.obs["total_guide_umis"] = np.asarray(counts.sum(axis=1)).ravel()
    guide_mod.obs["n_guides_detected"] = np.diff(counts.indptr)
    guide_mod.obs["n_guides_assigned"] = np.diff(calls.indptr)
    guide_mod.uns["guide_assignment"] = {
        "method": "pysceptre.assign_grnas_thresholding",
        "threshold": THRESHOLD,
    }

    gene.X = narrow(gene.X)
    for layer in list(gene.layers):
        gene.layers[layer] = narrow(gene.layers[layer])
    guide_mod.X = narrow(guide_mod.X)
    cut.update()

    out = OUT / "channels" / f"{ms}.h5mu"
    out.parent.mkdir(parents=True, exist_ok=True)
    partial = out.with_suffix(".partial.h5mu")
    cut.write_h5mu(partial, compression="gzip")
    partial.rename(out)

    knee_path = OUT / "knee_curves" / f"{ms}.tsv.gz"
    knee_path.parent.mkdir(parents=True, exist_ok=True)
    pd.concat(curves, ignore_index=True).to_csv(knee_path, sep="\t", index=False)

    held = np.asarray(cut.obsm["guide"]).ravel()
    return {
        "measurement_set": ms,
        "channel_name": info.entity_id,
        "sample": info.sample,
        "replicate": int(info.replicate),
        "day": int(info.day),
        "chip": info.chip,
        "batch": int(info.batch),
        "channel": int(info.channel),
        "barcodes_gene": curves[0]["n_barcodes"].iloc[0] if len(curves[0]) else 0,
        "barcodes_guide": curves[1]["n_barcodes"].iloc[0] if len(curves[1]) else 0,
        "cells_floor": int(cut.n_obs),
        "cells_with_guide_counts": int(held.sum()),
        "cells_assigned": int((np.diff(calls.indptr) > 0).sum()),
        "median_gene_umis": float(np.median(fractions["total_counts"])),
        "median_genes": float(np.median(fractions["n_genes_by_counts"])),
        "file_mb": round(out.stat().st_size / 1e6, 1),
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--channels", nargs="*")
    args = parser.parse_args(argv)

    t = table()
    if args.channels:
        t = t[t.measurement_sets.isin(args.channels)]
    summary_path = OUT / "channels.tsv"
    OUT.mkdir(parents=True, exist_ok=True)
    done = pd.read_csv(summary_path, sep="\t") if summary_path.exists() else pd.DataFrame()

    for info in t.itertuples():
        ms = info.measurement_sets
        if (OUT / "channels" / f"{ms}.h5mu").exists():
            continue
        start = time.time()
        row = prepare(info, args.workdir)
        row["seconds"] = round(time.time() - start)
        done = pd.concat(
            [done[done.get("measurement_set", pd.Series(dtype=str)) != ms], pd.DataFrame([row])],
            ignore_index=True,
        )
        done.to_csv(summary_path, sep="\t", index=False)
        shutil.rmtree(args.workdir / ms, ignore_errors=True)
        print(
            f"{info.sample} {ms} day {row['day']}: {row['cells_floor']:,} cells >= {FLOOR},"
            f" {row['cells_assigned']:,} with a guide, {row['file_mb']} MB, {row['seconds']} s",
            flush=True,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
