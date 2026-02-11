# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "marimo>=0.19.0",
#     "pyzmq>=27.1.0",
#     "scanpy>=1.12",
# ]
# ///

import marimo

__generated_with = "0.19.7"
app = marimo.App()


@app.cell
def _():
    import marimo as mo
    return


@app.cell
def _():
    import scanpy as sc
    import anndata as ad
    import pandas as pd
    return pd, sc


@app.cell
def _(pd, sc):
    def read_mtx_singlecol_features(path: str) -> sc.AnnData:
        adata = sc.read_mtx(f"{path}/matrix.mtx.gz").T  # cells x genes

        genes = pd.read_csv(f"{path}/features.tsv.gz", header=None, sep="\t")
        adata.var_names = genes[0].astype(str).values

        barcodes = pd.read_csv(f"{path}/barcodes.tsv.gz", header=None, sep="\t")
        adata.obs_names = barcodes[0].astype(str).values

        adata.var_names_make_unique()
        return adata


    return (read_mtx_singlecol_features,)


@app.cell
def _(read_mtx_singlecol_features, sc):
    samples = {
        "d0_ipsc": "../../data/external/10x_multiome_5_timepoints/rna/d0_ipsc/matrix",
        "d1_ps": "../../data/external/10x_multiome_5_timepoints/rna/d1_ps/matrix",
        "d2_meso": "../../data/external/10x_multiome_5_timepoints/rna/d2_meso/matrix",
        "d3_ec": "../../data/external/10x_multiome_5_timepoints/rna/d3_ec/matrix",
        "d4_ec": "../../data/external/10x_multiome_5_timepoints/rna/d4_ec/matrix",
    }

    adatas = {}
    for name, path in samples.items():
        a = read_mtx_singlecol_features(path)   # or read_10x_mtx(... var_names="gene_ids")
        a.obs["sample"] = name
        a.obs["timepoint"] = name.split("_")[0]  # d0, d1, ...
        a.obs["cell_type"] = name.split("_")[1]  # ipsc, ps, ...
        adatas[name] = a

    adata_all = sc.concat(adatas, label="sample", index_unique="-")
    adata_all.raw = None
    adata_all
    return (adata_all,)


@app.cell
def _(adata_all):
    adata_all.write_h5ad("../../results/2026-01-19_10x_multiome/rna/combined_5timepoints_cdna.h5ad")
    return


if __name__ == "__main__":
    app.run()
