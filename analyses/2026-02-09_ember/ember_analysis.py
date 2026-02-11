# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "ember-py>=0.1.1",
#     "marimo>=0.19.0",
#     "pyzmq>=27.1.0",
# ]
# ///

import marimo

__generated_with = "0.19.8"
app = marimo.App()


@app.cell
def _():
    import marimo as mo

    return


@app.cell
def _():
    # Import required packages
    import ember_py
    print(ember_py.__version__)
    from ember_py.light_ember import light_ember
    from ember_py.generate_pvals import generate_pvals
    from ember_py.plots import plot_partition_specificity, plot_block_specificity, plot_sample_counts, plot_psi_blocks
    from ember_py.top_genes import highly_specific_to_block, highly_specific_to_partition, non_specific_to_partition

    import scanpy as sc
    import gc
    import pandas as pd

    return generate_pvals, light_ember, plot_sample_counts, sc


@app.cell
def _(sc):
    adata = sc.read_h5ad('../../results/2026-01-19_10x_multiome/rna/combined_5timepoints_cdna.h5ad', backed = 'r')
    adata
    return (adata,)


@app.cell
def _(plot_sample_counts):
    plot_sample_counts(h5ad_dir = '../../results/2026-01-19_10x_multiome/rna/combined_5timepoints_cdna.h5ad',
        save_dir = '../../results/2026-02-09_ember',
        sample_id_col = 'sample',
        category_col = 'cell_type',
        condition_col = 'cell_type'
        )
    return


@app.cell
def _(adata):
    adata.obs['sample'].value_counts()
    return


@app.cell
def _(light_ember):
    light_ember(
        h5ad_dir = '../../results/2026-01-19_10x_multiome/rna/combined_5timepoints_cdna.h5ad',
        partition_label = 'timepoint',
        save_dir = '../../results/2026-02-09_ember',
        sampling=False,
        sample_id_col= 'sample',
        num_draws=100,
        partition_pvals=True,
        block_pvals=False,
        n_cpus=2
    )
    return


@app.cell
def _(generate_pvals):
    generate_pvals(h5ad_dir = '../../results/2026-01-19_10x_multiome/rna/combined_5timepoints_cdna.h5ad',
                   partition_label = 'timepoint',
                   entropy_metrics_dir = '../../results/2026-02-09_ember',
                   save_dir = '../../results/2026-02-09_ember',
                   sample_id_col = 'sample',
                   category_col = 'timepoint',
                   condition_col = 'timepoint',
                   n_iterations=100,
                   n_cpus=2)
    return


if __name__ == "__main__":
    app.run()
