# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "gcsfs>=2026.1.0",
#     "marimo>=0.19.0",
#     "pyzmq>=27.1.0",
#     "scanpy>=1.11.5",
# ]
# ///

import marimo

__generated_with = "0.19.4"
app = marimo.App()


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell
def _(mo):
    mo.md(r"""
    # Format the hashing metadata for CMO data

    **CMO** stands for **Cell Multiplexing Oligos**, which are used in 10x multiome experiments to label cells from different samples with unique oligonucleotide barcodes. This allows multiple samples to be pooled together in a single sequencing run, enabling efficient use of sequencing resources and reducing batch effects.

    This notebook takes the already processed CMO data present in oak and formats the metadata to be used in downstream analyses.

    Input folders (one for each 10x channel) from oak:
     - /oak/stanford/groups/engreitz/Users/dulguun/220519_EC_Profiling/flip_index/CMO1_CITE-seq-Count_outputs
     - /oak/stanford/groups/engreitz/Users/dulguun/220519_EC_Profiling/flip_index/CMO2_CITE-seq-Count_outputs

    Each of these folders have a subfolder for each sequencing `Lane` with the UMI counts in matrix market format(MTX). This is the general structure: `CMO{1,2}_CITE-seq-Count_outputs/Lane_{1,2,3,4}/umi_count`.
    Unfortunately, the feature file cannot be read in scanpy because of the formatting.

    The notebook will use the hashing metadata file present in `gs://fc-76565551-28c3-4b6d-a048-e272103bcbd1/metadata/10x_multiome_5_timepoints/cmo_metadata_complete.tsv` to fix the feature file, than read the MTX using scanpy and save the final AnnData object with formatted metadata.

    Output files:
     - One single h5ad data containing the counts for both 10x channels and the formatted metadata for CMO hashing.

     channel1 = `gs://fc-76565551-28c3-4b6d-a048-e272103bcbd1/data/10x_multiome_5_timepoints/cmo_channel1.annotated.h5ad`
     channel2 = `gs://fc-76565551-28c3-4b6d-a048-e272103bcbd1/data/10x_multiome_5_timepoints/cmo_channel2.annotated.h5ad`
    """)
    return


@app.cell
def _():
    import gzip
    import pandas as pd
    import numpy as np
    import scipy.sparse as sp
    import scanpy as sc
    from pathlib import Path
    return Path, gzip, np, pd, sc, sp


@app.cell
def _(mo):
    mo.md(r"""
    ### Reading the metadata file
    """)
    return


@app.cell
def _(pd):
    cmo_metadata_gcs = "gs://fc-76565551-28c3-4b6d-a048-e272103bcbd1/metadata/10x_multiome_5_timepoints/cmo_metadata_complete.tsv"

    cmo_metadata_df = pd.read_csv(cmo_metadata_gcs, sep="\t")

    cmo_metadata_df.head()
    return (cmo_metadata_df,)


@app.cell
def _(mo):
    mo.md(r"""
    ### Read the feature file and fix it

    I download the two main folders locally to read and fix the feature files.

    Change the paths below to the local paths where you downloaded the data or the oak paths if you want to read directly from oak.
    """)
    return


@app.cell
def _(Path, cmo_metadata_df, update_features_tsv):
    # ---- paths ----
    channel1_dir = Path("../../data/external/10x_multiome_5_timepoints/CMO1_CITE-seq-Count_outputs")
    channel2_dir = Path("../../data/external/10x_multiome_5_timepoints/CMO2_CITE-seq-Count_outputs")

    update_features_tsv(
        channel1_dir,
        cmo_metadata_df
    )
    update_features_tsv(
        channel2_dir,
        cmo_metadata_df
    )

    return channel1_dir, channel2_dir


@app.cell
def _(channel1_dir, sc):
    # Test reading in data
    test_adata = sc.read_10x_mtx(
        channel1_dir / "Lane_1/umi_count",
        var_names="gene_symbols",
        gex_only=False
    )

    test_adata.var.head(), test_adata.obs.head()
    return


@app.cell
def _(add_suffix, channel1_dir, sum_illumina_lanes):
    # Loading all lanes and combining for channel 1
    cmo_ch1 = sum_illumina_lanes(channel1_dir)
    cmo_ch1 = add_suffix(cmo_ch1, "_10x_5timepoints_channel1")
    print("channel1:", cmo_ch1.shape, "unique?", cmo_ch1.obs_names.is_unique)

    return (cmo_ch1,)


@app.cell
def _(cmo_ch1):
    print("channel1:", cmo_ch1.shape), cmo_ch1.obs.head()
    return


@app.cell
def _(add_suffix, channel2_dir, sum_illumina_lanes):
    # Loading all lanes and combining for channel 2
    cmo_ch2 = sum_illumina_lanes(channel2_dir)
    cmo_ch2 = add_suffix(cmo_ch2, "_10x_5timepoints_channel1")
    print("channel2:", cmo_ch2.shape, "unique?", cmo_ch2.obs_names.is_unique)
    return (cmo_ch2,)


@app.cell
def _(cmo_ch2):
    print("channel2:", cmo_ch2.shape), cmo_ch2.obs.head()
    return


@app.cell
def _(mo):
    mo.md(r"""
    ### Annotate the best CMO
    """)
    return


@app.cell
def _(Path, annotate_top_cmo, cmo_ch1, cmo_ch2, write_anndata):
    cmo_ch1_annot = annotate_top_cmo(cmo_ch1, prefix="cmo")
    cmo_ch2_annot = annotate_top_cmo(cmo_ch2, prefix="cmo")

    # Write the results

    out_dir = Path("../../data/external/10x_multiome_5_timepoints")
    write_anndata(cmo_ch1_annot, out_dir / "cmo_channel1.annotated.h5ad")
    write_anndata(cmo_ch2_annot, out_dir / "cmo_channel2.annotated.h5ad")
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## Helper functions
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ### Updating the feature files helpers
    """)
    return


@app.cell
def _(Path, gzip, pd):
    def update_feature_tsv(
        mtx_dir: Path,
        metadata: pd.DataFrame,
        feature_type: str = "Antibody Capture",
        make_backup: bool = True,
    ) -> pd.DataFrame:
        """
        Update features.tsv.gz in a 10x MTX folder to be 10x-compatible (3 columns),
        using Biosample as feature_name and Hashtag as feature_id.

        Expects original features.tsv.gz to be 1-column entries like:
            CMO1-TTGTCACGGTAATTA
        """
        features_path = mtx_dir / "features.tsv.gz"
        if not features_path.exists():
            raise FileNotFoundError(f"{features_path} not found")

        # read features header to check column count
        with gzip.open(features_path, "rt") as f:
            first_line = f.readline().rstrip("\n")
        n_cols = first_line.count("\t") + 1

        if n_cols >= 3:
            # Already 10x-style; do nothing
            print(f"Skipping {features_path} (already {n_cols} columns)")
            return None

        # build lookup from provided metadata df
        required = {"Barcode", "Hashtag", "Biosample"}
        missing = required - set(metadata.columns)
        if missing:
            raise ValueError(f"metadata is missing columns: {sorted(missing)}")

        lookup = metadata.loc[:, ["Barcode", "Hashtag", "Biosample"]].set_index("Barcode")

        # read features
        with gzip.open(features_path, "rt") as f:
            feats_raw = pd.read_csv(f, header=None, sep="\t")

        if feats_raw.shape[1] != 1:
            raise ValueError(
                f"Expected 1-column features.tsv.gz, found {feats_raw.shape[1]} columns in {features_path}"
            )

        feats_raw.columns = ["raw"]
        feats_raw["Barcode"] = feats_raw["raw"].astype(str).str.split("-").str[-1]

        # map to metadata
        mapped = feats_raw.join(lookup, on="Barcode")

        out = pd.DataFrame({
            0: mapped["Hashtag"].fillna(mapped["raw"]).astype(str),      # feature_id
            1: mapped["Biosample"].fillna(mapped["raw"]).astype(str),    # feature_name
            2: feature_type,                                             # feature_type
        })

        # backup
        if make_backup:
            backup_path = features_path.with_suffix(".bak")
            if not backup_path.exists():
                features_path.rename(backup_path)
            else:
                raise FileExistsError(f"Backup already exists: {backup_path}")

        # write updated features
        with gzip.open(features_path, "wt") as f:
            out.to_csv(f, sep="\t", header=False, index=False)

        return out


    def update_features_tsv(base_path: Path, metadata: pd.DataFrame, lane_glob: str = "Lane_*") -> None:
        """Update features.tsv.gz for all Lane_*/umi_count folders under base_path."""
        for lane_dir in sorted(base_path.glob(lane_glob)):
            lane_umi = lane_dir / "umi_count"
            if not lane_umi.exists():
                continue
            print(f"Updating features in {lane_umi} ...")
            update_feature_tsv(lane_umi, metadata)
    return (update_features_tsv,)


@app.cell
def _(mo):
    mo.md(r"""
    ### Combining the counts helpers
    """)
    return


@app.cell
def _(Path, np, pd, sc, sp):
    def load_cmo_lane_adata(lane_umi_dir: Path) -> sc.AnnData:
        """Load one Illumina lane umi_count folder as AnnData (keeps CMO features)."""
        return sc.read_10x_mtx(lane_umi_dir, var_names="gene_symbols", gex_only=False)

    def sum_illumina_lanes(base_dir: Path, lane_glob: str = "Lane_*") -> sc.AnnData:
        """
        Sum counts across Illumina lanes (Lane_1..Lane_4) inside one CITE-seq-Count output dir.
        Returns a single AnnData with union of barcodes and summed X.
        """
        lane_ads = []
        for lane_dir in sorted(base_dir.glob(lane_glob)):
            lane_umi = lane_dir / "umi_count"
            if lane_umi.exists():
                lane_ads.append(load_cmo_lane_adata(lane_umi))

        if not lane_ads:
            raise FileNotFoundError(f"No {lane_glob}/umi_count found under {base_dir}")

        # Ensure identical var order across lanes
        ref_vars = lane_ads[0].var_names
        for i, ad in enumerate(lane_ads):
            if not ad.var_names.equals(ref_vars):
                lane_ads[i] = ad[:, ref_vars].copy()

        # Union of barcodes
        all_barcodes = np.array(sorted(set().union(*[set(ad.obs_names) for ad in lane_ads])))
        bc_to_row = {bc: i for i, bc in enumerate(all_barcodes)}

        # Build one matrix per lane aligned to union barcodes, then sum (avoids CSR structure edits)
        mats = []
        for ad in lane_ads:
            row_idx = np.fromiter((bc_to_row[bc] for bc in ad.obs_names),
                                  dtype=np.int64, count=ad.n_obs)
            # Permutation/placement matrix: (union_cells x lane_cells)
            P = sp.csr_matrix(
                (np.ones(ad.n_obs, dtype=np.int8), (row_idx, np.arange(ad.n_obs))),
                shape=(len(all_barcodes), ad.n_obs),
            )
            mats.append(P @ ad.X)  # (union_cells x features)

        X_sum = mats[0]
        for M in mats[1:]:
            X_sum = X_sum + M

        out = sc.AnnData(X_sum, obs=pd.DataFrame(index=all_barcodes), var=lane_ads[0].var.copy())
        out.obs_names = all_barcodes
        out.var_names = ref_vars
        return out

    def add_suffix(adata: sc.AnnData, suffix: str) -> sc.AnnData:
        ad = adata.copy()
        ad.obs_names = [f"{bc}{suffix}" for bc in ad.obs_names.astype(str)]
        return ad
    return add_suffix, sum_illumina_lanes


@app.cell
def _(mo):
    mo.md(r"""
    ### Annotate best CMO helpers
    """)
    return


@app.cell
def _(Path, np, pd, sc):
    def annotate_top_cmo(
        adata: sc.AnnData,
        prefix: str = "cmo",
        store_second: bool = True,
    ) -> sc.AnnData:
        """
        Add per-cell top CMO annotations to adata.obs.

        Writes (with the given prefix):
          - {prefix}_best
          - {prefix}_best_count
          - {prefix}_second        (optional)
          - {prefix}_second_count  (optional)
          - {prefix}_margin        (optional: best - second)

        Returns a *copy* of the AnnData with annotations added.
        """
        ad = adata.copy()

        X = ad.X
        if hasattr(X, "toarray"):  # sparse
            X = X.toarray()

        if X.shape[1] == 0:
            raise ValueError("AnnData has 0 variables; cannot compute top CMO.")

        if not store_second:
            imax = X.argmax(axis=1)
            best = X[np.arange(X.shape[0]), imax]
            ad.obs[f"{prefix}_best"] = pd.Categorical(ad.var_names[imax].astype(str))
            ad.obs[f"{prefix}_best_count"] = best
            return ad

        # top-2 per row
        top2_idx = np.argpartition(X, -2, axis=1)[:, -2:]  # (n_cells, 2), unsorted
        top2_val = np.take_along_axis(X, top2_idx, axis=1)

        # sort within each row so col1 is best
        order = np.argsort(top2_val, axis=1)
        top2_idx = np.take_along_axis(top2_idx, order, axis=1)
        top2_val = np.take_along_axis(top2_val, order, axis=1)

        second_idx = top2_idx[:, 0]
        best_idx   = top2_idx[:, 1]
        second_val = top2_val[:, 0]
        best_val   = top2_val[:, 1]

        ad.obs[f"{prefix}_best"] = pd.Categorical(ad.var_names[best_idx].astype(str))
        ad.obs[f"{prefix}_best_count"] = best_val
        ad.obs[f"{prefix}_second"] = pd.Categorical(ad.var_names[second_idx].astype(str))
        ad.obs[f"{prefix}_second_count"] = second_val
        ad.obs[f"{prefix}_margin"] = best_val - second_val

        return ad


    def write_anndata(adata: sc.AnnData, out_path: Path) -> Path:
        """Write AnnData to .h5ad, creating parent dirs if needed."""
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(out_path)
        return out_path
    return annotate_top_cmo, write_anndata


if __name__ == "__main__":
    app.run()
