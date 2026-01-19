# /// script
# requires-python = ">=3.12"
# dependencies = [
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
    return


@app.cell
def _():
    import scanpy as sc
    return (sc,)


@app.cell
def _(sc):
    adata = sc.read_h5ad("../../data/external/10x_multiome_5_timepoints/10x_5timepoints_channel1.h5ad")

    adata

    return (adata,)


@app.cell
def _(adata):
    # mitochondrial genes, "MT-" for human, "Mt-" for mouse
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    # hemoglobin genes
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
    return


@app.cell
def _(adata, sc):
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)
    return


@app.cell
def _(adata, sc):
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
    )
    return


@app.cell
def _(adata, sc):
    sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

    return


@app.cell
def _(adata, sc):
    sc.pp.filter_cells(adata, min_counts=1000)

    return


@app.cell
def _(adata):
    adata.shape
    return


@app.cell
def _(adata, sc):
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
    )
    return


@app.cell
def _(cmo_data, sc):
    cmo1_fnp = "../../data/external/10x_multiome_5_timepoints/CMO1_CITE-seq-Count_outputs/Lane_1/umi_count"

    cmo_adata = sc.read_10x_mtx(cmo1_fnp)
    cmo_data
    return


@app.cell
def _(adata):
    adata.obs.index
    return


@app.cell
def _():
    return


@app.cell
def _():
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
