# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "altair>=6.0.0",
#     "marimo>=0.19.0",
#     "pyzmq>=27.1.0",
#     "snapatac2>=2.8.0",
#     "umap-learn>=0.5.11",
# ]
# ///

import marimo

__generated_with = "0.19.7"
app = marimo.App()


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell
def _():
    import altair as alt
    import numpy as np
    import snapatac2 as snap
    import glob
    import umap
    from pathlib import Path
    return Path, alt, glob, np, snap


@app.cell
def _(Path, glob):
    # Get the list of fragment files
    paths = sorted(map(Path, glob.glob("../../data/external/10x_multiome_5_timepoints/atac/d*.fragments.bed.gz")))

    # Create a list of tuples (timepoint_name, file_path)
    files = [(p.name.split("_atac")[0], p) for p in paths]
    # e.g. [('d0_ipsc', Path('...d0_ipsc_atac.sorted.fragments.bed.gz')), ...]
    files

    outpath = "../../results/2026-01-19_10x_multiome/atac/"
    return files, outpath


@app.cell
def _(files, outpath, snap):
    adatas = snap.pp.import_fragments(
        [fl for _, fl in files],
        file=[outpath + name + '.h5ad' for name, _ in files],
        chrom_sizes=snap.genome.hg38,
        sorted_by_barcode=False,
    )
    return (adatas,)


@app.cell
def _(adatas, snap):
    snap.metrics.tsse(adatas, snap.genome.hg38)
    snap.pp.filter_cells(adatas, min_tsse=7)
    snap.pp.add_tile_matrix(adatas, bin_size=1000)
    snap.pp.select_features(adatas, n_features=50000)
    return


@app.cell
def _(adatas):
    adatas
    return


@app.cell
def _(adatas, files, outpath, snap):
    data = snap.AnnDataSet(
        adatas=[(name, adata) for (name, _), adata in zip(files, adatas)],
        filename=f"{outpath}combined_5timepoints_atac.h5ad"
    )
    data
    return (data,)


@app.cell
def _(data, snap):
    peak_file_fnp = "../../data/external/10x_multiome_5_timepoints/scE2G/all_timepoints_merged_regions-macs2_peaks.narrowPeak.sorted.candidateRegions.bed"
    peak_matrix = snap.pp.make_peak_matrix(
        data,
        peak_file=peak_file_fnp
    )
    peak_matrix
    return (peak_matrix,)


@app.cell
def _(outpath, peak_matrix):
    peak_matrix.write(f'{outpath}/atac_peak_count_matrix.h5ad')
    return


@app.cell
def _(data, np):
    print(f'Number of cells: {data.n_obs}')
    print(f'Number of unique barcodes: {np.unique(data.obs_names).size}')
    return


@app.cell
def _(data, np):
    #unique_cell_ids = [sa + ':' + bc for sa, bc in zip(data.obs['sample'], data.obs_names)]
    #data.obs_names = unique_cell_ids
    assert data.n_obs == np.unique(data.obs_names).size
    return


@app.cell
def _(data, snap):
    snap.pp.select_features(data, n_features=50000)
    return


@app.cell
def _(data, snap):
    snap.tl.spectral(data)
    return


@app.cell
def _(data, snap):
    snap.tl.umap(data)
    return


@app.cell
def _(data, mo, plot_embedding_altair):
    chart = plot_embedding_altair(data, color="sample")
    mo.ui.altair_chart(chart)
    return


@app.cell
def _(data, snap):
    snap.pp.knn(data)
    snap.tl.leiden(data)
    return


@app.cell
def _(data, mo, plot_embedding_altair):
    _chart = plot_embedding_altair(data, color="leiden")
    mo.ui.altair_chart(_chart)
    return


@app.cell
def _(data, mo, plot_cluster_sample_proportions):
    _chart = plot_cluster_sample_proportions(data, cluster="leiden", sample="sample")
    mo.ui.altair_chart(_chart)
    return


@app.cell
def _(data, mo, plot_cluster_sample_heatmap):
    _chart = plot_cluster_sample_heatmap(
        data,
        cluster="leiden",
        sample="sample",
        color_scheme="greens",
        cell_size=24
    )
    mo.ui.altair_chart(_chart)
    return


@app.cell
def _(data):
    data.close()
    return


@app.cell
def _(alt, pd):
    def plot_cluster_sample_heatmap(
        adata,
        cluster="leiden",
        sample="sample",
        sort_clusters=True,
        sort_samples=True,
        color_scheme="reds",
        color_domain=None,   # e.g. (0, 0.5) or [0, 0.8]
        width=500,
        cell_size=22,
        show_text=True,
    ):
        # Pull obs columns robustly into pandas
        df = pd.DataFrame(index=adata.obs_names)
        for colname in [cluster, sample]:
            col = adata.obs[colname]
            try:
                df[colname] = col.to_numpy()
            except Exception:
                try:
                    df[colname] = col.to_list()
                except Exception:
                    df[colname] = list(col)

        # counts and proportions within cluster
        tab = df.groupby([cluster, sample]).size().reset_index(name="n")
        tab["prop"] = tab["n"] / tab.groupby(cluster)["n"].transform("sum")

        # Ordering
        if sort_clusters:
            try:
                cluster_order = sorted(df[cluster].astype(int).unique())
                cluster_order = [str(x) for x in cluster_order]
                tab[cluster] = tab[cluster].astype(int).astype(str)
            except Exception:
                cluster_order = sorted(df[cluster].astype(str).unique())
                tab[cluster] = tab[cluster].astype(str)
        else:
            cluster_order = None
            tab[cluster] = tab[cluster].astype(str)

        if sort_samples:
            sample_order = df[sample].astype(str).value_counts().index.tolist()
            tab[sample] = tab[sample].astype(str)
        else:
            sample_order = None
            tab[sample] = tab[sample].astype(str)

        # Build scale safely (only include domain if provided)
        scale_kwargs = {"scheme": color_scheme}
        if color_domain is not None:
            scale_kwargs["domain"] = list(color_domain)
            scale_kwargs["clamp"] = True

        heat = (
            alt.Chart(tab)
            .mark_rect()
            .encode(
                x=alt.X(
                    f"{sample}:N",
                    title="Sample",
                    sort=sample_order,
                    axis=alt.Axis(labelAngle=-45),
                ),
                y=alt.Y(
                    f"{cluster}:N",
                    title="Cluster",
                    sort=cluster_order,
                ),
                color=alt.Color(
                    "prop:Q",
                    title="Proportion",
                    scale=alt.Scale(**scale_kwargs),
                    legend=alt.Legend(format=".0%"),
                ),
                tooltip=[
                    alt.Tooltip(f"{cluster}:N", title="Cluster"),
                    alt.Tooltip(f"{sample}:N", title="Sample"),
                    alt.Tooltip("n:Q", title="Cells"),
                    alt.Tooltip("prop:Q", title="Proportion", format=".1%"),
                ],
            )
            .properties(width=alt.Step(cell_size), height=alt.Step(cell_size))
        )

        if not show_text:
            return heat.properties(width=width)

        text = (
            alt.Chart(tab)
            .mark_text()
            .encode(
                x=alt.X(f"{sample}:N", sort=sample_order),
                y=alt.Y(f"{cluster}:N", sort=cluster_order),
                text=alt.Text("prop:Q", format=".0%"),
            )
        )

        return (heat + text).properties(width=width)
    return (plot_cluster_sample_heatmap,)


@app.cell
def _():
    return


@app.cell
def _(alt, pd):
    def plot_cluster_sample_proportions(
        adata,
        cluster="leiden",
        sample="sample",
        sort_clusters=True,
        width=700,
        height=300,
    ):
        # Pull obs columns robustly into pandas
        df = pd.DataFrame(index=adata.obs_names)
        for colname in [cluster, sample]:
            col = adata.obs[colname]
            try:
                df[colname] = col.to_numpy()
            except Exception:
                try:
                    df[colname] = col.to_list()
                except Exception:
                    df[colname] = list(col)

        # counts and proportions within cluster
        counts = (
            df.groupby([cluster, sample])
              .size()
              .reset_index(name="n")
        )
        counts["prop"] = counts["n"] / counts.groupby(cluster)["n"].transform("sum")

        # optional: nicer cluster ordering (numeric if possible)
        if sort_clusters:
            try:
                # if clusters are like "0","1","2",... this makes them sort numerically
                counts[cluster] = counts[cluster].astype(int).astype(str)
            except Exception:
                pass

        chart = (
            alt.Chart(counts)
            .mark_bar()
            .encode(
                x=alt.X(f"{cluster}:N", title="Cluster"),
                y=alt.Y("prop:Q", title="Proportion", axis=alt.Axis(format="%")),
                color=alt.Color(f"{sample}:N", title="Sample"),
                tooltip=[
                    alt.Tooltip(f"{cluster}:N", title="Cluster"),
                    alt.Tooltip(f"{sample}:N", title="Sample"),
                    alt.Tooltip("n:Q", title="Cells"),
                    alt.Tooltip("prop:Q", title="Proportion", format=".1%"),
                ],
            )
            .properties(width=width, height=height)
        )

        return chart
    return (plot_cluster_sample_proportions,)


@app.cell
def _(alt):
    import pandas as pd

    def plot_embedding_altair(
        adata,
        basis="X_umap",          # e.g. "X_umap", "X_spectral"
        color="sample",          # adata.obs column
        dims=(0, 1),             # which columns of the embedding to plot
        point_size=8,
        opacity=0.6,
        width=600,
        height=600,
        max_points=None,
        seed=0,
    ):
        emb = adata.obsm[basis]
        i, j = dims

        df = pd.DataFrame(
            {f"{basis}_{i+1}": emb[:, i], f"{basis}_{j+1}": emb[:, j]},
            index=adata.obs_names,
        )

        col = adata.obs[color]
        # robust conversion for pandas / polars / arrow-backed objects
        try:
            df[color] = col.to_numpy()
        except Exception:
            try:
                df[color] = col.to_list()
            except Exception:
                df[color] = list(col)

        if max_points is not None and len(df) > max_points:
            df = df.sample(n=max_points, random_state=seed)

        df = df.reset_index(names="cell")

        xcol = f"{basis}_{i+1}"
        ycol = f"{basis}_{j+1}"

        return (
            alt.Chart(df)
            .mark_circle(size=point_size, opacity=opacity)
            .encode(
                x=alt.X(f"{xcol}:Q", title=f"{basis} {i+1}"),
                y=alt.Y(f"{ycol}:Q", title=f"{basis} {j+1}"),
                color=alt.Color(f"{color}:N", legend=alt.Legend(title=color)),
                tooltip=["cell:N", alt.Tooltip(f"{color}:N"), alt.Tooltip(f"{xcol}:Q"), alt.Tooltip(f"{ycol}:Q")],
            )
            .properties(width=width, height=height)
        )
    return pd, plot_embedding_altair


if __name__ == "__main__":
    app.run()
