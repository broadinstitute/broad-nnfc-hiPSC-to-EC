# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "anndata==0.13.4",
#     "crispr-analysis-utils",
#     "h5py>=3.12.0",
#     "jupyter-scatter==1.0.1",
#     "marimo>=0.25.2",
#     "matplotlib==3.11.2",
#     "mudata==0.4.1",
#     "numpy==2.5.5",
#     "pandas==3.1.1",
#     "pysceptre[fast,io]==0.6.0",
#     "scipy==1.18.1",
# ]
#
# [tool.uv.sources]
# crispr-analysis-utils = { git = "https://github.com/EngreitzLab/crispr-analysis-utils", rev = "6e75a186e32ebcf479b8fc5d4440ea8610848c7d" }
# ///

import marimo

__generated_with = "0.25.2"
app = marimo.App(width="medium", auto_download=["html"])


@app.cell(hide_code=True)
def title(mo):
    mo.md(r"""
    # Day-2 RNA and guide QC

    Day 2 of the hiPSC-to-endothelial-cell differentiation screen, read out on a
    targeted panel of 292 genes, from its day object. The steps are the pysceptre
    tutorial's (`analyses/2026-11-17_pysceptre_tutorial`), through guide QC,
    guide assignment, the pairs and the covariates; the calibration check, power
    and discovery come later, in their own notebook.
    """)
    return


@app.cell
def import_marimo():
    import marimo as mo

    return (mo,)


@app.cell
def imports(mo):
    from pathlib import Path

    import anndata as ad
    import h5py
    import matplotlib
    import matplotlib.pyplot as plt
    import mudata
    import numpy as np
    import pandas as pd
    import pysceptre
    from anndata.io import read_elem, sparse_dataset
    from crispr_analysis_utils.guide_annotation.db import load_store
    from crispr_analysis_utils.pairs import gene_tss
    from scipy import sparse, stats

    # marimo shows figures as SVG, with labels as real text: sharp at any zoom
    plt.rcParams["savefig.format"] = "svg"
    plt.rcParams["svg.fonttype"] = "none"
    plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial", "Helvetica"]

    mo.md(f"pysceptre **{pysceptre.__version__}**, mudata **{mudata.__version__}**")
    return (
        Path,
        ad,
        gene_tss,
        h5py,
        load_store,
        matplotlib,
        mudata,
        np,
        pd,
        plt,
        pysceptre,
        read_elem,
        sparse,
        sparse_dataset,
        stats,
    )


@app.cell
def jscatter_import():
    import jscatter

    return (jscatter,)


@app.cell
def settings(mo):
    # Cells to analyse, a random sample so a training run finishes in minutes.
    # None analyses every cell; only then do the numbers match the full screen.
    N_CELLS = None
    SAMPLE_SEED = 0

    mo.md(
        f"Analysing a random **{N_CELLS:,}** cells: set `N_CELLS = None` for all of them."
        if N_CELLS is not None
        else "Analysing every cell."
    )
    return N_CELLS, SAMPLE_SEED


@app.cell(hide_code=True)
def md_get_data(mo):
    mo.md(r"""
    ## Get the data

    The day-2 object, `results/2026-11-12_day_objects/days/day2.h5mu`, holds every
    barcode past its channel's 100-UMI floor, with the QC flags `flag_cells.py`
    set in `obs`. The GENCODE v43 store the guide design was annotated with gives
    each gene its TSS.
    """)
    return


@app.cell(hide_code=True)
def get_data(Path, mo):
    data_path = mo.notebook_dir().parents[1] / "results" / "2026-11-12_day_objects" / "days" / "day2.h5mu"
    STORE = Path.home() / "Annotations" / "HG38" / "GENCODE" / "gencode.v43.annotation.parquet"
    mo.md(
        f"- `{data_path}`: {data_path.stat().st_size / 1e9:.2f} GB\n"
        f"- `{STORE}`: {STORE.stat().st_size / 1e6:.0f} MB"
    )
    return STORE, data_path


@app.cell(hide_code=True)
def md_filter(mo):
    mo.md(r"""
    ## Keep the cells that pass QC

    `flag_cells.py` flagged every barcode and dropped none. `pass_rna`: past the
    day's gene UMI cut and inside each channel's MAD bands for gene UMIs, genes
    detected and mitochondrial share. `pass_guide`: at least one guide assigned at
    20 UMIs. `pass_qc`: both.

    On top of `pass_qc` we apply the three day-2 limits of the R QC
    (`analyses/2026-01-20_load_counts_and_qc`): at least 1,000 gene UMIs, at least
    45% of them from readout genes (`pct_counts_targeted`) and at least 16% from
    ribosomal genes (`pct_counts_ribo_rmd`, the R QC's genes). The day object
    records the last two in `uns["qc_rules"]` but does not apply them. We keep the
    cells that pass all of it (`keep`).
    """)
    return


@app.cell
def qc_cells(data_path, h5py, mo, pd, read_elem):
    # the R QC's limits (analyses/2026-01-20_load_counts_and_qc), applied on top of pass_qc
    MIN_GENE_UMIS, MIN_READOUT_PCT, MIN_RIBO_PCT = 1_000, 45, 16

    with h5py.File(data_path, "r") as _f:
        obs = read_elem(_f["obs"])
        _flags = read_elem(_f["uns/qc_flags"])
        _gene_obs = read_elem(_f["mod/gene/obs"])

    _metrics = _gene_obs[["total_counts", "pct_counts_targeted", "pct_counts_ribo_rmd"]].reindex(obs.index)
    obs["pass_limits"] = (  # a barcode missing from the gene modality compares False
        (_metrics["total_counts"] >= MIN_GENE_UMIS)
        & (_metrics["pct_counts_targeted"] >= MIN_READOUT_PCT)
        & (_metrics["pct_counts_ribo_rmd"] >= MIN_RIBO_PCT)
    )
    obs["keep"] = obs["pass_qc"] & obs["pass_limits"]

    # the pass_qc cells' metrics, to show what the limits drop
    cell_metrics = _metrics[obs["pass_qc"].to_numpy(dtype=bool)]
    assert cell_metrics.notna().all().all(), "a pass_qc cell has no gene metrics"

    mo.vstack(
        [
            mo.md(
                f"Gene UMI cut **{_flags['gene_umi_cut']:,}**, MAD bands at **{_flags['nmads']}** MADs "
                f"on {', '.join(_flags['band'])}. Limits: **{MIN_GENE_UMIS:,}** gene UMIs, "
                f"**{MIN_READOUT_PCT}%** readout, **{MIN_RIBO_PCT}%** ribosomal."
            ),
            pd.Series(
                {
                    "barcodes": len(obs),
                    "past the gene UMI cut": int(obs["rna_umi_cut"].sum()),
                    "pass_rna": int(obs["pass_rna"].sum()),
                    "pass_guide": int(obs["pass_guide"].sum()),
                    "pass_qc": int(obs["pass_qc"].sum()),
                    f"pass_qc, under {MIN_GENE_UMIS:,} gene UMIs": int((cell_metrics["total_counts"] < MIN_GENE_UMIS).sum()),
                    f"pass_qc, under {MIN_READOUT_PCT}% readout": int((cell_metrics["pct_counts_targeted"] < MIN_READOUT_PCT).sum()),
                    f"pass_qc, under {MIN_RIBO_PCT}% ribosomal": int((cell_metrics["pct_counts_ribo_rmd"] < MIN_RIBO_PCT).sum()),
                    "keep: pass_qc and all three limits": int(obs["keep"].sum()),
                },
                name="cells",
            ).to_frame(),
        ]
    )
    return MIN_GENE_UMIS, MIN_READOUT_PCT, MIN_RIBO_PCT, cell_metrics, obs


@app.cell
def gene_tss_table(STORE, gene_tss, load_store, mo):
    tss = gene_tss(load_store(STORE))
    mo.md(f"**{len(tss):,} genes** with a TSS in `{STORE.name}`")
    return (tss,)


@app.cell(hide_code=True)
def md_load(mo):
    mo.md(r"""
    ## Load the screen

    Read only the cells kept, straight from the file, so the barcodes that
    fail are never loaded. Two modalities over the same cells: `gene` (UMIs per
    gene) and `guide` (UMIs per guide, plus a `guide_assignment` layer saying
    which guides each cell carries). The day object's three kb layers stay behind.

    Each gene gains `gene_tss` (0-based) and `gene_strand`: the start of its
    MANE Select transcript, else its Ensembl canonical one, else its longest CDS.
    """)
    return


@app.cell
def load(
    N_CELLS,
    SAMPLE_SEED,
    ad,
    data_path,
    h5py,
    mo,
    mudata,
    np,
    obs,
    pd,
    read_elem,
    sparse_dataset,
    tss,
):
    _rows = np.flatnonzero(obs["keep"].to_numpy(dtype=bool))
    _n_all = _rows.size
    if N_CELLS is not None and N_CELLS < _n_all:
        # sorted, so the sampled cells keep their order in the file; the same cells every run
        _rows = _rows[np.sort(np.random.default_rng(SAMPLE_SEED).choice(_n_all, N_CELLS, replace=False))]
    _cells = obs.index[_rows]

    def _modality(f, name, layers=()):
        """One modality, restricted to the kept cells and to the layers named."""
        group = f["mod"][name]
        mod_obs = read_elem(group["obs"])
        rows = np.flatnonzero(mod_obs.index.isin(_cells))
        if rows.size != len(_cells):
            raise ValueError(f"{name}: {len(_cells) - rows.size:,} kept cells are not in this modality")
        return ad.AnnData(
            X=sparse_dataset(group["X"])[rows],
            obs=mod_obs.iloc[rows],
            var=read_elem(group["var"]),
            layers={layer: sparse_dataset(group["layers"][layer])[rows] for layer in layers},
            uns=read_elem(group["uns"]),
        )

    with h5py.File(data_path, "r") as _f:
        _gene = _modality(_f, "gene")
        _guide = _modality(_f, "guide", layers=("guide_assignment",))
        _uns = read_elem(_f["uns"])

    # gene var_names carry the Ensembl version; the store's ids do not
    _ids = _gene.var_names.str.replace(r"\.\d+$", "", regex=True)
    _gene.var["gene_tss"] = pd.array(tss["tss"].reindex(_ids), dtype="Int64")
    _gene.var["gene_strand"] = pd.Categorical(tss["strand"].reindex(_ids))
    _no_tss = _gene.var["gene_tss"].isna() & _gene.var["readout"]
    if _no_tss.any():
        raise ValueError(f"{int(_no_tss.sum())} readout genes have no TSS in the GENCODE store")

    mdata = mudata.MuData({"gene": _gene, "guide": _guide})
    mdata.obs = obs.reindex(mdata.obs_names)
    mdata.uns.update(_uns)

    mo.vstack(
        [
            mo.md(
                f"**{mdata.n_obs:,}** of the {_n_all:,} cells kept"
                + (", a random sample (`N_CELLS`)" if mdata.n_obs < _n_all else "")
                + f". {int(_gene.var['gene_tss'].notna().sum()):,} of {_gene.n_vars:,} genes have a TSS, "
                "every readout gene among them."
            ),
            mo.plain_text(repr(mdata)),
        ]
    )
    return (mdata,)


@app.cell(hide_code=True)
def md_guides(mo):
    mo.md(r"""
    ## The guide library

    `mdata["guide"].var` is the per-guide table, in the IGVF format. For each
    guide it says what kind of guide it is (`type`), which element it targets
    (`intended_target_name` and its coordinates), and which genes the design
    expects it to act on (`putative_target_genes`), where it names any.

    An element is identified by its coordinates, not its name: a promoter is named
    by its gene (`ENSG…`) and an exon by its exon (`ENSE…`), so several elements
    can share one name.
    """)
    return


@app.cell
def guides(mdata, mo):
    guides = mdata["guide"].var
    _elements = (
        guides[["intended_target_chr", "intended_target_start", "intended_target_end"]]
        .dropna()
        .drop_duplicates()
    )

    mo.vstack(
        [
            mo.md(
                f"**{len(guides):,} guides** against **{len(_elements):,} elements**"
            ),
            guides["type"].value_counts().rename("guides").to_frame(),
            mo.ui.table(
                guides[
                    [
                        "spacer",
                        "type",
                        "genomic_element",
                        "intended_target_name",
                        "intended_target_chr",
                        "intended_target_start",
                        "intended_target_end",
                        "putative_target_genes",
                    ]
                ].reset_index(names="guide_id")
            ),
        ]
    )
    return (guides,)


@app.cell(hide_code=True)
def md_genes(mo):
    mo.md(r"""
    ## The genes

    `mdata["gene"].var` describes every gene of the GENCODE v43 annotation the
    counts were made against, and `mdata.uns["reference"]` records which file that
    was. Only the `readout` genes are captured by the
    targeted panel, so those are the genes we can test. `gene_tss` and
    `gene_strand` give one TSS per gene, which we need to pair elements with
    nearby genes.
    """)
    return


@app.cell
def genes(mdata, mo, pd):
    genes = mdata["gene"].var

    mo.vstack(
        [
            mo.md(
                f"**{len(genes):,} genes**, **{int(genes['readout'].sum()):,}** of them "
                "in the readout panel"
            ),
            pd.Series(mdata.uns["reference"], name="value").to_frame(),
            mo.ui.table(
                genes.loc[
                    genes["readout"],
                    [
                        "symbol",
                        "gene_type",
                        "gene_chr",
                        "gene_start",
                        "gene_end",
                        "gene_tss",
                        "gene_strand",
                    ],
                ].reset_index(names="gene_id")
            ),
        ]
    )
    return (genes,)


@app.cell(hide_code=True)
def md_rna_qc(mo):
    mo.md(r"""
    ## RNA QC: the knee plot

    Rank each channel's barcodes by their gene UMIs. Cells sit on the high
    plateau; empty droplets, holding only ambient RNA, form the low shoulder
    beneath it. The cut the day object applied (`uns["qc_flags"]["gene_umi_cut"]`)
    should fall on the cliff between the two.

    We kept only the cells that passed QC, so the whole curves come from
    `uns["knee_curves"]`: every channel's barcodes, sampled at log-spaced ranks.
    """)
    return


@app.cell
def rna_knee(mdata, plot_knees):
    _curves = mdata.uns["knee_curves"]
    plot_knees(
        _curves[_curves["modality"] == "gene"],
        cut=mdata.uns["qc_flags"]["gene_umi_cut"],
        ylabel="Gene UMIs per barcode",
        title=f"The {mdata.uns['qc_flags']['gene_umi_cut']:,}-UMI cut falls on every channel's cliff",
        subtitle="Day 2: barcodes ranked by gene UMIs, log scales; dashed line, the cut",
    )
    return


@app.cell(hide_code=True)
def md_readout(mo):
    mo.md(r"""
    ## Do the UMIs come from the readout genes?

    DC-TAP-seq amplifies the 292 readout genes with gene-specific primers, so most
    of each cell's gene UMIs should come from them. We compute every cell's
    readout share from the counts (readout UMIs over all gene UMIs) and check it
    against the `pct_counts_targeted` column the QC stored.
    """)
    return


@app.cell
def readout_share(genes, mdata, mo, np, pd, plot_share):
    _counts = mdata["gene"].X
    _readout_umis = np.asarray(_counts @ genes["readout"].to_numpy(dtype=float)).ravel()
    _total_umis = np.asarray(_counts.sum(axis=1)).ravel()
    readout_share = pd.Series(
        _readout_umis / _total_umis, index=mdata["gene"].obs_names, name="readout_share"
    )
    _agrees = np.allclose(100 * readout_share, mdata["gene"].obs["pct_counts_targeted"])

    mo.vstack(
        [
            mo.md(
                f"Pooled over all cells, the readout genes hold "
                f"**{_readout_umis.sum() / _total_umis.sum():.1%}** of the gene UMIs. "
                f"Matches the stored `pct_counts_targeted`: **{_agrees}**."
            ),
            plot_share(
                readout_share,
                mark=0.5,
                xlabel="Share of the cell's gene UMIs from readout genes, %",
                title=f"Readout genes supply {readout_share.median():.0%} of a typical cell's UMIs",
                subtitle=(
                    f"Day 2, {len(readout_share):,} cells: share of each cell's gene UMIs "
                    f"from the {int(genes['readout'].sum())} readout genes"
                ),
            ),
        ]
    )
    return (readout_share,)


@app.cell(hide_code=True)
def md_readout_by_batch(mo):
    mo.md(r"""
    ### By batch

    The same share, one row per batch, each row as a percentage of its own cells,
    so the small batch is not flattened by the large ones.
    """)
    return


@app.cell
def readout_share_by_batch(genes, mdata, plot_share_by, readout_share):
    _batch = mdata.obs.loc[readout_share.index, "batch"]
    _medians = readout_share.groupby(_batch.to_numpy()).median()
    plot_share_by(
        readout_share,
        _batch,
        xlim=(50, 100),
        xlabel="Share of the cell's gene UMIs from readout genes, %",
        title=f"Batch medians run from {_medians.min():.0%} to {_medians.max():.0%}",
        subtitle=(
            f"Day 2: share of each cell's gene UMIs from the {int(genes['readout'].sum())} "
            f"readout genes; {int((readout_share < 0.5).sum())} cells below 50% not shown"
        ),
    )
    return


@app.cell(hide_code=True)
def md_readout_by_channel(mo):
    mo.md(r"""
    ### By channel

    Pick a batch to see its channels, one row each.
    """)
    return


@app.cell
def batch_picker(mdata, mo):
    batch_picker = mo.ui.dropdown(
        options=sorted(map(str, mdata.obs["batch"].unique())), value="4", label="Batch"
    )
    batch_picker
    return (batch_picker,)


@app.cell
def readout_share_by_channel(
    batch_picker,
    mdata,
    plot_share_by,
    readout_share,
):
    _obs = mdata.obs.loc[readout_share.index]
    _in_batch = (_obs["batch"] == batch_picker.value).to_numpy()
    _channels = _obs.loc[_in_batch, "channel_name"]
    _medians = readout_share[_in_batch].groupby(_channels.to_numpy()).median()
    plot_share_by(
        readout_share[_in_batch],
        _channels,
        label="",
        row_height=0.6,
        xlim=(50, 100),
        xlabel="Share of the cell's gene UMIs from readout genes, %",
        title=(
            f"Batch {batch_picker.value}: channel medians run from "
            f"{_medians.min():.0%} to {_medians.max():.0%}"
        ),
        subtitle="Day 2: share of each cell's gene UMIs from the readout genes, one row per channel",
    )
    return


@app.cell(hide_code=True)
def md_umi_vs_ribo(mo):
    mo.md(r"""
    ### Gene UMIs against the ribosomal share

    Each cell's gene UMIs (`total_counts`) against the share of them from ribosomal
    genes (`pct_counts_ribo`), both as the day object stores them. The plot is
    interactive: scroll to zoom, drag to pan, and hover over a point for its values.
    """)
    return


@app.cell
def umi_vs_ribo(mdata, pd, threshold_scatter):
    _obs = mdata["gene"].obs[["total_counts", "pct_counts_ribo"]]
    _tenth = pd.qcut(_obs["total_counts"], 10, labels=False)
    _fewest = _obs.loc[_tenth == 0, "pct_counts_ribo"].median()
    _most = _obs.loc[_tenth == 9, "pct_counts_ribo"].median()
    # the median of each twentieth of the cells, ranked by gene UMIs
    _medians = _obs.groupby(pd.qcut(_obs["total_counts"], 20, labels=False)).median()
    threshold_scatter(
        _obs,
        "total_counts",
        "pct_counts_ribo",
        xlabel="Gene UMIs per cell",
        ylabel="Ribosomal share of gene UMIs, %",
        line=zip(_medians["total_counts"].tolist(), _medians["pct_counts_ribo"].tolist()),
        title=f"The ribosomal share {'falls' if _most < _fewest else 'rises'} from {_fewest:.0f}% to {_most:.0f}% as gene UMIs rise",
        subtitle="Day 2: each cell's gene UMIs and the share of them from ribosomal genes, one point per cell",
        caption=(
            f"Red line: the median ribosomal share in each twentieth of the {len(_obs):,} cells kept, ranked by gene UMIs. "
            f"The median share is {_fewest:.1f}% in the tenth of cells with the fewest UMIs and "
            f"{_most:.1f}% in the tenth with the most; Spearman correlation "
            f"{_obs['total_counts'].corr(_obs['pct_counts_ribo'], method='spearman'):.2f}."
        ),
    )
    return


@app.cell(hide_code=True)
def md_mito_ribo(mo):
    mo.md(r"""
    ### Mitochondrial and ribosomal shares

    The R QC (`analyses/2026-01-20_load_counts_and_qc`) drew these plots for every
    barcode before it filtered; here they are for the cells kept. The ribosomal
    share is `pct_counts_ribo_rmd`, the R QC's own definition: the genes matching
    `^RPS|^RPL|^RRP|^MRP|^FAU`. Its third histogram, the share of gene UMIs from
    readout genes, is the readout-share plot above.
    """)
    return


@app.cell
def mito_ribo_hist(mdata, mo, plot_hist_median):
    _obs = mdata["gene"].obs
    _mito, _ribo = _obs["pct_counts_mito"], _obs["pct_counts_ribo_rmd"]
    mo.vstack(
        [
            plot_hist_median(
                _mito,
                xlabel="Mitochondrial share of gene UMIs, %",
                title=f"A median cell has {_mito.median():.2f}% of its UMIs from mitochondrial genes",
                subtitle="Day 2: each cell's mitochondrial share, in 50 bins",
            ),
            plot_hist_median(
                _ribo,
                xlabel="Ribosomal share of gene UMIs, %",
                title=f"A median cell has {_ribo.median():.1f}% of its UMIs from ribosomal genes",
                subtitle="Day 2: each cell's ribosomal share, the R QC's genes, in 50 bins",
            ),
            mo.md(f"Bars: number of cells, of the {len(_obs):,} cells kept. Line: the median."),
        ]
    )
    return


@app.cell(hide_code=True)
def md_umi_vs_readout(mo):
    mo.md(r"""
    ### What the R QC's limits drop

    The `pass_qc` cells, before the three limits applied when the cells were chosen.
    The first plot draws the UMI and readout limits: a cell below either line is
    dropped. The second shows only the cells below the ribosomal limit, which are
    dropped whatever their UMIs and readout share.

    Both plots are interactive: scroll to zoom, drag to pan, and hover over a point
    for its values. The marginal histograms of the R plots are not drawn.
    """)
    return


@app.cell
def umi_vs_readout(
    MIN_GENE_UMIS,
    MIN_READOUT_PCT,
    cell_metrics,
    threshold_scatter,
):
    _few_umis = cell_metrics["total_counts"] < MIN_GENE_UMIS
    _low_readout = cell_metrics["pct_counts_targeted"] < MIN_READOUT_PCT
    threshold_scatter(
        cell_metrics,
        "total_counts",
        "pct_counts_targeted",
        MIN_GENE_UMIS,
        MIN_READOUT_PCT,
        xlabel="Gene UMIs per cell",
        ylabel="Readout share of gene UMIs, %",
        title=(
            f"The UMI and readout limits drop {int((_few_umis | _low_readout).sum()):,} "
            f"of the {len(cell_metrics):,} pass_qc cells"
        ),
        subtitle="Day 2: each pass_qc cell's gene UMIs and the share of them from readout genes, one point per cell",
        caption=(
            f"Lines: the limits, {MIN_GENE_UMIS:,} gene UMIs and {MIN_READOUT_PCT}% from readout genes. "
            f"{int(_few_umis.sum()):,} cells have fewer UMIs and {int(_low_readout.sum()):,} a lower readout share."
        ),
    )
    return


@app.cell
def umi_vs_readout_low_ribo(
    MIN_GENE_UMIS,
    MIN_READOUT_PCT,
    MIN_RIBO_PCT,
    cell_metrics,
    mo,
    threshold_scatter,
):
    _low = cell_metrics[cell_metrics["pct_counts_ribo_rmd"] < MIN_RIBO_PCT]
    _only_ribo = (_low["total_counts"] >= MIN_GENE_UMIS) & (_low["pct_counts_targeted"] >= MIN_READOUT_PCT)
    (
        threshold_scatter(
            _low,
            "total_counts",
            "pct_counts_targeted",
            MIN_GENE_UMIS,
            MIN_READOUT_PCT,
            xlabel="Gene UMIs per cell",
            ylabel="Readout share of gene UMIs, %",
            title=f"The ribosomal limit drops {len(_low):,} pass_qc cells, {int(_only_ribo.sum()):,} of them inside the other two",
            subtitle=f"Day 2: the pass_qc cells under {MIN_RIBO_PCT}% ribosomal UMIs, gene UMIs against the readout share",
            caption=f"Lines: the limits, {MIN_GENE_UMIS:,} gene UMIs and {MIN_READOUT_PCT}% from readout genes.",
        )
        if len(_low)
        else mo.md(f"No pass_qc cell has under {MIN_RIBO_PCT}% ribosomal UMIs.")
    )
    return


@app.cell(hide_code=True)
def md_gene_upm(mo):
    mo.md(r"""
    ### Readout genes against the rest

    For every gene, two means over the cells kept: UMIs per cell, and UPM, the UMIs
    per million of the cell's gene UMIs. As in the R QC, the readout genes are
    compared with all the others on UPM, on a log scale, with 1 added so that a
    gene never seen stays on the axis.
    """)
    return


@app.cell
def gene_means(mdata, mo, np, pd, plot_upm_box):
    # per gene, over the cells kept: mean UMIs per cell, and mean UMIs per million of the cell's gene UMIs
    _X = mdata["gene"].X
    _totals = np.asarray(_X.sum(axis=1)).ravel()
    gene_means = pd.DataFrame(
        {
            "symbol": mdata["gene"].var["symbol"].to_numpy(),
            "readout": mdata["gene"].var["readout"].to_numpy(),
            "mean_umi": np.asarray(_X.mean(axis=0)).ravel(),
            "mean_upm": (_X.T @ (1e6 / _totals)) / _X.shape[0],
        },
        index=mdata["gene"].var_names,
    )
    _readout = gene_means[gene_means["readout"]]
    _other = gene_means[~gene_means["readout"]]
    mo.vstack(
        [
            plot_upm_box(
                gene_means["mean_upm"],
                gene_means["readout"],
                title=(
                    f"Readout genes have a median {_readout['mean_upm'].median():,.0f} UPM, "
                    f"the rest {_other['mean_upm'].median():.3f}"
                ),
                subtitle="Day 2: each gene's mean UPM over the cells kept, plus 1, log scale",
            ),
            mo.md(
                f"Boxes: the middle half of the genes, the line their median; whiskers reach 1.5 times "
                f"the box, and points beyond are single genes. {len(_readout):,} readout genes, "
                f"{len(_other):,} others, {int((_other['mean_upm'] == 0).sum()):,} of them never seen."
            ),
        ]
    )
    return (gene_means,)


@app.cell(hide_code=True)
def md_markers(mo):
    mo.md(r"""
    ### Marker and housekeeping genes

    The panel flags some readout genes as markers of day 0, 2 or 4, or as
    housekeeping genes (`metadata/metadata_genes_final.csv`). A gene with more than
    one flag counts once, in that order, as in the R QC.
    """)
    return


@app.cell
def marker_upm(gene_means, mo, np, pd, plot_marker_strip):
    _meta = pd.read_csv(mo.notebook_dir().parents[1] / "metadata" / "metadata_genes_final.csv")
    _meta.index = _meta["intended_target_gene_id"].str.replace(r"\.\d+$", "", regex=True)
    _meta = _meta[~_meta.index.duplicated()]  # a gene with two primers has two rows
    _flags = {
        "is_d0_marker": "Day 0 marker",
        "is_d2_marker": "Day 2 marker",
        "is_d4_marker": "Day 4 marker",
        "is_housekeeping_gene": "Housekeeping",
    }
    _flagged = _meta[list(_flags)].astype(str).apply(lambda _s: _s.str.lower().eq("true"))
    _class = pd.Series(np.select([_flagged[_f] for _f in _flags], list(_flags.values()), ""), index=_meta.index)
    _class = _class[_class != ""]

    _upm = pd.Series(gene_means["mean_upm"].to_numpy(), index=gene_means.index.str.replace(r"\.\d+$", "", regex=True))
    assert _class.index.isin(_upm.index).all(), "a flagged gene is not in the gene modality"
    _upm = _upm.reindex(_class.index)
    _medians = _upm.groupby(_class).median()
    mo.vstack(
        [
            plot_marker_strip(
                _upm,
                _class,
                list(_flags.values()),
                title=(
                    f"Day-0 markers have a median {_medians['Day 0 marker']:,.0f} UPM; "
                    f"day-2 {_medians['Day 2 marker']:,.1f}, day-4 {_medians['Day 4 marker']:,.1f}"
                ),
                subtitle="Day 2: each flagged gene's mean UPM over the cells kept, log scale",
            ),
            mo.md(
                "Points: one per gene, spread sideways so they don't overlap. Red line and number: "
                f"the class's median. {int((_upm == 0).sum())} flagged genes have 0 UPM and are not drawn."
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def md_reference(mo):
    mo.md(r"""
    ### Against a reference day-0 dataset

    Each readout gene's mean UMIs per cell here, against its UPM in the day-0 10x
    5′ pilot in `data/external/WTC11_10x_5prime`, joined on the Ensembl id, with
    the R QC's lines at 30 UPM and 1 UMI. A gene the reference never saw is drawn
    at 0.01 UPM. The R QC also set that gene's UMIs here to 0.01, which put it in
    Q3 whatever it measured; here the quadrants are counted on the values
    themselves. Hover over a point for its gene.
    """)
    return


@app.cell
def reference_comparison(gene_means, mo, pd, threshold_scatter):
    REF_UPM_CUT, MEAN_UMI_CUT = 30, 1  # the R QC's lines
    _ref = pd.read_csv(
        mo.notebook_dir().parents[1]
        / "data"
        / "external"
        / "WTC11_10x_5prime"
        / "D0_10x5prime_pilot_counts_summary_ribo15_mt15_umi2000.tsv",
        sep="\t",
        index_col="ensembl_id",
    )
    assert _ref.index.is_unique
    _readout = gene_means.loc[gene_means["readout"], ["symbol", "mean_umi"]]
    _readout.index = _readout.index.str.replace(r"\.\d+$", "", regex=True)
    _both = _readout.join(_ref["mean_upm"].rename("reference_upm"), how="inner")
    _high_ref, _high_here = _both["reference_upm"] >= REF_UPM_CUT, _both["mean_umi"] >= MEAN_UMI_CUT
    _quadrants = pd.DataFrame(
        {
            "reference_upm": [f">= {REF_UPM_CUT}", f"< {REF_UPM_CUT}", f"< {REF_UPM_CUT}", f">= {REF_UPM_CUT}"],
            "mean_umi_here": [f">= {MEAN_UMI_CUT}", f">= {MEAN_UMI_CUT}", f"< {MEAN_UMI_CUT}", f"< {MEAN_UMI_CUT}"],
            "genes": [
                int((_high_ref & _high_here).sum()),
                int((~_high_ref & _high_here).sum()),
                int((~_high_ref & ~_high_here).sum()),
                int((_high_ref & ~_high_here).sum()),
            ],
        },
        index=pd.Index(["Q1", "Q2", "Q3", "Q4"], name="quadrant"),
    )
    _shown = _both.assign(reference_upm=_both["reference_upm"].clip(lower=0.01))
    mo.vstack(
        [
            threshold_scatter(
                _shown,
                "reference_upm",
                "mean_umi",
                REF_UPM_CUT,
                MEAN_UMI_CUT,
                xlabel="UPM in the reference",
                ylabel="Mean UMIs per cell here",
                scale="log",
                tooltip=("x", "y", "symbol"),
                title=(
                    "Readout genes' UMIs here follow their reference UPM "
                    f"(Spearman {_both['reference_upm'].corr(_both['mean_umi'], method='spearman'):.2f})"
                ),
                subtitle="Day 2: each readout gene's mean UMIs per cell here and UPM in the reference, log scales",
                caption=(
                    f"{len(_both)} of {len(_readout)} readout genes are in the reference; "
                    f"{int((_both['reference_upm'] < 0.01).sum())} it never saw are drawn at 0.01 UPM."
                ),
            ),
            _quadrants,
        ]
    )
    return


@app.cell(hide_code=True)
def md_guide_qc(mo):
    mo.md(r"""
    ## Guide QC: the knee plot

    The same ranking, for guide UMIs. If the same droplets hold cells, the guide
    knee falls at the same rank as the gene knee.
    """)
    return


@app.cell
def guide_knee(mdata, plot_knees):
    _curves = mdata.uns["knee_curves"]
    plot_knees(
        _curves[_curves["modality"] == "guide"],
        ylabel="Guide UMIs per barcode",
        title="Guide UMIs fall off the same cliff as gene UMIs",
        subtitle="Day 2: barcodes ranked by guide UMIs, log scales",
    )
    return


@app.cell(hide_code=True)
def md_assignment(mo):
    mo.md(r"""
    ## Assign guides to cells

    A cell carries a guide when it holds at least a threshold of that guide's
    UMIs. The day object called guides at 20 UMIs and kept only the cells with at
    least one, so this file can show only what a *higher* threshold changes: how
    many cells keep a guide, and how many guides each carries (the MOI), per
    replicate. Pick a threshold with the slider; pysceptre's thresholding method
    then makes the calls.
    """)
    return


@app.cell
def threshold_sweep(guides_at, mdata, pd):
    _counts = mdata["guide"].X
    _replicate = mdata.obs.loc[mdata["guide"].obs_names, "replicate"].to_numpy()
    _rows = []
    for _threshold in [20, 30, 40, 50, 75, 100, 150, 200]:
        _by = pd.Series(guides_at(_counts, _threshold)).groupby(_replicate)
        _rows.append(
            pd.DataFrame(
                {
                    "threshold": _threshold,
                    "replicate": _by.size().index,
                    "cells": _by.size().to_numpy(),
                    "with_guide": _by.apply(lambda n: (n > 0).mean()).to_numpy(),
                    "moi": _by.mean().to_numpy(),
                }
            )
        )
    threshold_sweep = pd.concat(_rows, ignore_index=True)
    threshold_sweep
    return (threshold_sweep,)


@app.cell
def guide_threshold(mo):
    guide_threshold = mo.ui.slider(
        start=20, stop=200, step=10, value=20, label="Guide UMI threshold", show_value=True
    )
    guide_threshold
    return (guide_threshold,)


@app.cell
def threshold_sweep_plot(guide_threshold, plot_sweep, threshold_sweep):
    plot_sweep(
        threshold_sweep,
        chosen=guide_threshold.value,
        title="A higher threshold lowers the MOI more than it loses cells",
        subtitle="Day 2: cells keeping a guide, and guides per cell, by guide UMI threshold",
    )
    return


@app.cell
def guide_assignment(guide_threshold, mdata, mo, np, pd, pysceptre):
    _guide = mdata["guide"]
    assignment = pysceptre.assign_grnas_thresholding(
        _guide.X.T, _guide.var_names, threshold=guide_threshold.value
    )
    _per_cell = pd.Series(np.asarray(assignment.assigned.sum(axis=0)).ravel())
    _by = _per_cell.groupby(mdata.obs.loc[_guide.obs_names, "replicate"].to_numpy())
    _summary = pd.DataFrame(
        {
            "cells": _by.size(),
            "with ≥ 1 guide, %": (100 * _by.apply(lambda n: (n > 0).mean())).round(1),
            "guides per cell, mean": _by.mean().round(2),
        }
    ).rename_axis("replicate")

    _check = ""
    if assignment.threshold == _guide.uns["guide_assignment"]["threshold"]:
        _differ = (assignment.assigned != _guide.layers["guide_assignment"].T.astype(bool)).nnz
        _check = f" Entries that differ from the stored layer: **{_differ}**."
    mo.vstack(
        [
            mo.md(
                f"Threshold **{assignment.threshold:.0f} UMIs**: "
                f"**{assignment.assigned.nnz:,}** guide-cell assignments.{_check}"
            ),
            _summary,
        ]
    )
    return (assignment,)


@app.cell(hide_code=True)
def md_guides_per_cell(mo):
    mo.md(r"""
    ## Guides per cell

    How many guides a cell carries is the screen's multiplicity of infection
    (MOI). With several guides per cell, a guide's controls are all the cells
    without it, which is pysceptre's `moi="high"` setting.
    """)
    return


@app.cell
def guides_per_cell(assignment, mdata, np, pd, plot_counts):
    guides_per_cell = pd.Series(
        np.asarray(assignment.assigned.sum(axis=0)).ravel(),
        index=mdata["guide"].obs_names,
        name="guides_per_cell",
    )
    plot_counts(
        guides_per_cell,
        xlabel="Guides assigned to the cell",
        title=f"The median cell carries {guides_per_cell.median():.0f} guides",
        subtitle=(
            f"Day 2, {len(guides_per_cell):,} cells: guides assigned per cell "
            f"({assignment.threshold:.0f}-UMI threshold); mean {guides_per_cell.mean():.1f}"
        ),
    )
    return (guides_per_cell,)


@app.cell(hide_code=True)
def note_low_moi(mo):
    mo.md(r"""
    /// note | Later: an internal low-MOI control

    About 30% of the cells carry exactly one called guide. Analysed on their
    own (`moi="low"`, against the cells whose only guide is non-targeting or
    safe-targeting), they are a low-MOI screen inside the high-MOI one; pairs
    called by both are robust. Caveats: one *called* guide is not one guide
    integrated, and weak guide capture is over-represented among these cells.
    ///
    """)
    return


@app.cell(hide_code=True)
def md_guides_vs_depth(mo):
    mo.md(r"""
    ## Guides per cell against depth

    A **doublet** is two cells captured in one droplet and read as one barcode. It
    carries both cells' guides, so it looks like a cell with many guides, and both
    cells' RNA, so it has about **twice the gene UMIs** of one cell. If the cells
    with many guides were mostly doublets, their gene UMIs would approach twice
    those of one-guide cells: the dashed line.

    These cells already passed QC, which removed cells more than five MADs from
    their channel's median gene UMIs, so the most obvious doublets may be gone: this
    checks what is left.
    """)
    return


@app.cell
def guides_vs_depth(assignment, guides_per_cell, mdata, np, pd, plot_depth):
    _cap = 16
    _genes = mdata["gene"].obs.loc[guides_per_cell.index, "total_counts"].to_numpy()
    _median = pd.Series(_genes).groupby(np.minimum(guides_per_cell.to_numpy(), _cap)).median()
    plot_depth(
        guides_per_cell.to_numpy(),
        _genes,
        cap=_cap,
        title=(
            "Cells with many guides are not twice as deep"
            if (_median / _median.loc[1]).max() < 1.5
            else "Cells with many guides are close to twice as deep"
        ),
        subtitle=(
            "Day 2: median gene UMIs per cell, relative to cells with one guide "
            f"({assignment.threshold:.0f}-UMI threshold)"
        ),
    )
    return


@app.cell(hide_code=True)
def md_threshold_calls(mo):
    mo.md(r"""
    ## Where the threshold sits

    Each guide detected in a cell has a UMI count, and most counts are background:
    ambient guide RNA, not a guide the cell carries. A guide call is a count at or
    above the threshold. Past the background, a well-placed threshold sits where
    moving it changes few calls.
    """)
    return


@app.cell
def threshold_calls(
    assignment,
    calls_by_threshold,
    mdata,
    mo,
    plot_calls_by_threshold,
):
    _threshold = int(assignment.threshold)
    _umis = mdata["guide"].X.data  # every nonzero guide count in every cell
    _calls = calls_by_threshold(_umis, top=max(200, 2 * _threshold))
    _lo, _hi = _threshold // 2, 2 * _threshold

    mo.vstack(
        [
            plot_calls_by_threshold(
                _calls,
                _threshold,
                title=f"From {_lo} to {_hi} UMIs, the threshold changes the calls by {1 - _calls[_hi] / _calls[_lo]:.0%}",
                subtitle="Day 2: guide calls, a guide at or above the threshold in a cell; shaded, half to twice it",
            ),
            mo.md(
                f"The curve starts at 2 UMIs: {(_umis == 1).mean():.0%} of the {len(_umis):,} nonzero guide "
                f"counts are a single UMI. A threshold of 1 would call {_calls[1] / mdata['guide'].n_obs:,.0f} "
                "guides per cell."
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def md_cells_per_guide(mo):
    mo.md(r"""
    ## Cells per guide

    A guide's cells are the group its effect is measured in, so a guide in few
    cells has little power. Count the cells assigned to each guide at the threshold
    chosen above, all 35 channels pooled.
    """)
    return


@app.cell
def cells_per_guide(assignment, mdata, np, pd, plot_cells_per):
    cells_per_guide = pd.Series(
        np.asarray(assignment.assigned.sum(axis=1)).ravel(),
        index=mdata["guide"].var_names,
        name="cells_assigned",
    )
    _none = int((cells_per_guide == 0).sum())
    plot_cells_per(
        cells_per_guide,
        unit="guide",
        bin_width=2,
        title=f"The median guide reaches {cells_per_guide.median():,.0f} cells",
        subtitle=(
            f"Day 2: cells assigned to each of the {len(cells_per_guide):,} guides "
            f"({assignment.threshold:.0f}-UMI threshold); {_none:,} guide{'s' * (_none != 1)} in no cell"
        ),
    )
    return (cells_per_guide,)


@app.cell(hide_code=True)
def md_cells_per_guide_by_type(mo):
    mo.md(r"""
    ### Cells per guide, by guide type

    Do the four guide types reach cells equally? The non-targeting guides feed the
    calibration check later, so they should be captured as well as the rest.
    """)
    return


@app.cell
def cells_per_guide_by_type(
    assignment,
    cells_per_guide,
    guides,
    plot_by_group,
):
    _types = guides["type"].astype(str).to_numpy()
    _median = cells_per_guide.groupby(_types).median()
    plot_by_group(
        cells_per_guide,
        _types,
        bin_width=2,
        unit="guides",
        xlabel="Cells assigned to the guide",
        title=f"{_median.idxmin().capitalize()} guides reach the fewest cells",
        subtitle=(
            f"Day 2: cells per guide by guide type, % of the type's guides "
            f"({assignment.threshold:.0f}-UMI threshold)"
        ),
    )
    return


@app.cell(hide_code=True)
def md_replicate_agreement(mo):
    mo.md(r"""
    ### Replicate agreement

    The replicates draw on the same guide library, so a guide with many UMIs in one
    should have many in the other. Each dot is a guide: its total UMIs in replicate
    1 against 2, and in 3 against 4, on log axes. A tight cloud means the
    replicates agree.
    """)
    return


@app.cell
def replicate_agreement(mdata, mo, np, pd, plot_replicates, sparse):
    _replicate = mdata.obs.loc[mdata["guide"].obs_names, "replicate"].astype(str).to_numpy()
    _names = sorted(set(_replicate))
    # cells × replicates: 1 where the cell belongs to the replicate
    _membership = sparse.csr_matrix(
        (np.ones(len(_replicate)), (np.arange(len(_replicate)), pd.Categorical(_replicate, categories=_names).codes)),
        shape=(len(_replicate), len(_names)),
    )
    # guides × replicates: every UMI of the guide in the replicate's cells
    guide_umis_by_replicate = pd.DataFrame(
        (mdata["guide"].X.T @ _membership).toarray(),
        index=mdata["guide"].var_names,
        columns=_names,
    )
    _pairs = [("1", "2"), ("3", "4")]
    _logs = np.log10(guide_umis_by_replicate.where(guide_umis_by_replicate > 0))
    _r = [_logs[_a].corr(_logs[_b]) for _a, _b in _pairs]

    mo.vstack(
        [
            plot_replicates(
                guide_umis_by_replicate,
                _pairs,
                title=(
                    f"Replicates {' and '.join(_pairs[int(np.argmax(_r))])} agree better than "
                    f"{' and '.join(_pairs[int(np.argmin(_r))])}"
                ),
                subtitle="Day 2: total UMIs of each guide in each replicate, one dot per guide",
            ),
            mo.md(
                "Median UMIs per guide, and guides with no UMI, per replicate: "
                + "; ".join(
                    f"**{_n}**: {guide_umis_by_replicate[_n].median():,.0f}, "
                    f"{int((guide_umis_by_replicate[_n] == 0).sum()):,} with none"
                    for _n in _names
                )
                + "."
            ),
        ]
    )
    return (guide_umis_by_replicate,)


@app.cell(hide_code=True)
def md_replicate_disagreement(mo):
    mo.md(r"""
    ### Where replicates 1 and 2 disagree

    Two places the disagreement could come from. **A channel**: compare each
    channel's guide UMIs with replicates 3 and 4 pooled, which agree with each
    other; a channel whose guides follow another profile stands out with a low
    correlation. **A set of guides**: take the guides whose UMIs differ more than
    10-fold between replicates 1 and 2, each scaled to its replicate's median, and
    count them by library (the guide name without its number) and guide type.
    """)
    return


@app.cell
def replicate_disagreement(guide_umis_by_replicate, mdata, mo, np, pd, sparse):
    _guide = mdata["guide"]
    _channel = mdata.obs.loc[_guide.obs_names, "channel_name"].astype(str).to_numpy()
    _names = sorted(set(_channel))
    # cells × channels: 1 where the cell belongs to the channel
    _membership = sparse.csr_matrix(
        (np.ones(len(_channel)), (np.arange(len(_channel)), pd.Categorical(_channel, categories=_names).codes)),
        shape=(len(_channel), len(_names)),
    )
    _by_channel = pd.DataFrame((_guide.X.T @ _membership).toarray(), index=_guide.var_names, columns=_names)
    _reference = np.log10(guide_umis_by_replicate["3"] + guide_umis_by_replicate["4"])


    def _r(umis):
        _both = (umis > 0) & np.isfinite(_reference)
        return np.corrcoef(np.log10(umis[_both]), _reference[_both])[0, 1]


    _obs = mdata.obs.drop_duplicates("channel_name").set_index("channel_name")
    channel_agreement = (
        pd.DataFrame(
            {
                "replicate": _obs.loc[_names, "replicate"].astype(str).to_numpy(),
                "cells": pd.Series(_channel).value_counts().reindex(_names).to_numpy(),
                "guides_with_umis": (_by_channel > 0).sum().to_numpy(),
                "r_vs_reps_3_4": [_r(_by_channel[_n]) for _n in _names],
            },
            index=pd.Index(_names, name="channel_name"),
        )
        .sort_values("r_vs_reps_3_4")
        .round({"r_vs_reps_3_4": 2})
    )

    # guides more than 10-fold apart between replicates 1 and 2, each replicate scaled to its median
    _scaled = {_n: (guide_umis_by_replicate[_n] + 1) / guide_umis_by_replicate[_n].median() for _n in ("1", "2")}
    _log_ratio = np.log10(_scaled["1"] / _scaled["2"])
    off_diagonal_guides = (
        pd.DataFrame(
            {
                "library": _guide.var_names.str.replace(r"_\d+$", "", regex=True),
                "type": _guide.var["type"].astype(str).to_numpy(),
                "higher_in_1": (_log_ratio > 1).to_numpy(),
                "higher_in_2": (_log_ratio < -1).to_numpy(),
            },
            index=_guide.var_names,
        )
        .groupby(["library", "type"])
        .agg(guides=("higher_in_1", "size"), higher_in_1=("higher_in_1", "sum"), higher_in_2=("higher_in_2", "sum"))
    )

    mo.vstack(
        [
            mo.md(
                f"**{int((_log_ratio > 1).sum()):,}** guides are more than 10-fold higher in replicate 1 and "
                f"**{int((_log_ratio < -1).sum()):,}** in replicate 2, of {len(_log_ratio):,}."
            ),
            mo.ui.table(channel_agreement),
            off_diagonal_guides,
        ]
    )
    return


@app.cell(hide_code=True)
def md_replicate_cells_per_guide(mo):
    mo.md(r"""
    ### Few cells per guide?

    Replicates 1 and 2 have fewer cells than 3 and 4, so each guide is carried by
    fewer cells there, and a guide's UMIs follow how many cells happened to carry
    it. If that is the cause, the replicates should agree once each guide is in
    enough cells in both. Count the cells each guide is assigned to in each
    replicate (the calls above), and compare replicates only over the guides with at
    least *k* cells in both.
    """)
    return


@app.cell
def replicate_cells_per_guide(
    assignment,
    guide_umis_by_replicate,
    mdata,
    mo,
    np,
    pd,
    sparse,
):
    _guide = mdata["guide"]
    _replicate = mdata.obs.loc[_guide.obs_names, "replicate"].astype(str).to_numpy()
    _names = sorted(set(_replicate))
    _membership = sparse.csr_matrix(
        (np.ones(len(_replicate)), (np.arange(len(_replicate)), pd.Categorical(_replicate, categories=_names).codes)),
        shape=(len(_replicate), len(_names)),
    )
    # guides × replicates: cells each guide is assigned to (assignment.assigned is guides × cells)
    _cells_per_guide = pd.DataFrame(
        (assignment.assigned.astype(float) @ _membership).toarray(), index=_guide.var_names, columns=_names
    )


    def _r(a, b, k):
        _both = (_cells_per_guide[a] >= k) & (_cells_per_guide[b] >= k)
        _x, _y = guide_umis_by_replicate.loc[_both, a], guide_umis_by_replicate.loc[_both, b]
        return int(_both.sum()), np.corrcoef(np.log10(_x), np.log10(_y))[0, 1]


    _rows = []
    for _k in (1, 3, 5, 10):
        for _a, _b in (("1", "2"), ("3", "4")):
            _n, _rk = _r(_a, _b, _k)
            _rows.append({"at least k cells in both": _k, "replicates": f"{_a} vs {_b}", "guides": _n, "r": round(_rk, 2)})

    mo.vstack(
        [
            mo.md(
                "Median cells per guide: "
                + ", ".join(f"replicate {_n} **{_cells_per_guide[_n].median():.0f}**" for _n in _names)
                + "."
            ),
            pd.DataFrame(_rows).set_index(["at least k cells in both", "replicates"]),
        ]
    )
    return


@app.cell
def md_guides_together(mo):
    mo.md(r"""
    ### Guides found together

    If a cell's guides arrive independently, two guides share a cell now and then by
    chance. How often depends on how common each guide is in the replicate and how
    many guides that replicate's cells carry. For every pair of guides seen in a cell
    together, we compare the cells they share with that chance count, and keep the
    pairs far above it: a Poisson test, Bonferroni-corrected over every possible pair.
    """)
    return


@app.cell
def guides_together(assignment, guides, mdata, mo, np, pd, sparse, stats):
    _calls = assignment.assigned.astype(np.int32)  # guides × cells, 1 where the guide is assigned
    _calls.eliminate_zeros()
    _k = np.asarray(_calls.sum(axis=0)).ravel()  # guides per cell
    _replicate = mdata.obs.loc[mdata["guide"].obs_names, "replicate"].astype(str).to_numpy()
    _names, _codes = np.unique(_replicate, return_inverse=True)
    _membership = sparse.csr_matrix(
        (np.ones(len(_codes)), (np.arange(len(_codes)), _codes)), shape=(len(_codes), len(_names))
    )

    # chance: each guide's share of its replicate's calls, times the ordered guide pairs in that replicate's cells
    _per_replicate = (_calls @ _membership).toarray()
    _share = _per_replicate / _per_replicate.sum(axis=0)
    _pairs_in_cells = np.bincount(_codes, weights=_k * (_k - 1.0))

    # every pair of guides seen together, and the cells they share
    _shared = sparse.triu(_calls @ _calls.T, k=1).tocoo()
    _expected = (_share[_shared.row] * _share[_shared.col] * _pairs_in_cells).sum(axis=1)
    _p = stats.poisson.sf(_shared.data - 1, _expected)
    _tests = _calls.shape[0] * (_calls.shape[0] - 1) / 2
    _keep = _p < 0.05 / _tests

    guides_together = pd.DataFrame(
        {
            "guide_a": guides.index[_shared.row[_keep]],
            "guide_b": guides.index[_shared.col[_keep]],
            "type_a": guides["type"].to_numpy()[_shared.row[_keep]],
            "type_b": guides["type"].to_numpy()[_shared.col[_keep]],
            "shared_cells": _shared.data[_keep],
            "expected": _expected[_keep],
            "p_value": _p[_keep],
        }
    ).sort_values("shared_cells", ascending=False, ignore_index=True)
    _a = guides.index.get_indexer(guides_together["guide_a"])
    _b = guides.index.get_indexer(guides_together["guide_b"])

    _element = (
        guides["intended_target_chr"].astype(str)
        + ":" + guides["intended_target_start"].astype(str)
        + "-" + guides["intended_target_end"].astype(str)
    ).where(guides["intended_target_chr"].notna()).to_numpy()
    _same_element = int(sum(isinstance(_x, str) and _x == _y for _x, _y in zip(_element[_a], _element[_b])))

    # spacers that differ only in their first base, and in their shared cells, the weaker one's UMIs over the stronger's
    _last19 = guides["spacer"].astype(str).str[-19:].to_numpy()
    _twin = _last19[_a] == _last19[_b]
    _library_twins = int(pd.Series(_last19).value_counts().pipe(lambda _c: _c * (_c - 1) // 2).sum())
    _spill = []
    for _i, _j in zip(_a[_twin], _b[_twin]):
        _umis = mdata["guide"].X[:, [_i, _j]].toarray()
        _both = _calls[[_i, _j]].toarray().all(axis=0)
        _spill.extend(_umis[_both].min(axis=1) / _umis[_both].max(axis=1))

    mo.vstack(
        [
            mo.md(
                f"**{len(guides_together):,}** pairs of guides share more cells than chance, involving "
                f"{len(set(_a) | set(_b)):,} of the {len(guides):,} guides. A typical pair shares "
                f"{guides_together['shared_cells'].median():.0f} cells where {guides_together['expected'].median():.2f} "
                f"are expected. Only {_same_element:,} of them target the same element, so sharing an element "
                "doesn't explain them.\n\n"
                f"{int(_twin.sum())} pairs are spacers that differ only in their first base (the library has "
                f"{_library_twins}). In their shared cells the weaker guide carries a median {np.median(_spill):.1%} "
                "of the stronger's UMIs: one guide's reads spill into its twin's count, and a cell with thousands "
                f"of UMIs passes the {assignment.threshold:.0f}-UMI threshold for both."
            ),
            mo.ui.table(guides_together, selection=None),
        ]
    )
    return (guides_together,)


@app.cell
def md_identical_guide_sets(mo):
    mo.md(r"""
    ### Cells with identical guide sets

    Where do the other pairs come from? Look at cells carrying exactly the same three
    or more guides. Chance can't make those: with tens of thousands of guides, two
    cells drawing the same set essentially never happens. Two explanations remain:

    - **Clones.** One transduced cell divided before the cells were loaded. A
      replicate's cells were split across its channels, so a clone's cells spread
      over the channels in proportion to their sizes.
    - **One cell, several barcodes.** A droplet with two or more beads reads one cell
      as several barcodes. These stay in one channel, and each barcode holds only
      part of the cell's RNA.

    For each group of cells with the same set, we ask whether all its cells sit in
    one channel, and how often channel sizes alone would put them there.

    Why it matters: copies of one cell are not independent observations, and guides
    that always travel together can't be told apart.
    """)
    return


@app.cell
def identical_guide_sets(
    assignment,
    guides,
    guides_together,
    mdata,
    mo,
    np,
    pd,
    plot_one_channel,
):
    _cells = assignment.assigned.T.tocsr()  # cells × guides
    _cells.eliminate_zeros()
    _cells.sort_indices()
    _k = np.diff(_cells.indptr)  # guides per cell
    _rows = np.flatnonzero(_k >= 3)
    _obs = mdata.obs.loc[mdata["guide"].obs_names, ["replicate", "channel_name"]].astype(str)

    # a cell's guide set as one key: its sorted guide positions, as bytes
    _sets = _obs.iloc[_rows].assign(
        key=[_cells.indices[_cells.indptr[_i] : _cells.indptr[_i + 1]].tobytes() for _i in _rows],
        barcode=lambda d: d.index.str.split("_").str[0],
    )
    _sets = _sets[_sets.groupby("key")["key"].transform("size") > 1]
    _groups = _sets.groupby("key", sort=False)

    # chance that a group's cells all sit in one channel: the replicate's channel shares, to the power of its cells
    _share = _obs.iloc[_rows].groupby("replicate")["channel_name"].value_counts(normalize=True)
    identical_sets = pd.DataFrame(
        {
            "cells": _groups.size(),
            "channels": _groups["channel_name"].nunique(),
            "replicate": _groups["replicate"].first(),
        }
    )
    identical_sets["one_channel_by_chance"] = [
        float((_share.loc[_r].to_numpy() ** _n).sum())
        for _r, _n in zip(identical_sets["replicate"], identical_sets["cells"])
    ]
    _members = [np.frombuffer(_key, dtype=_cells.indices.dtype) for _key in identical_sets.index]
    identical_sets["guides"] = [", ".join(guides.index[_m]) for _m in _members]
    # the guide pairs found together that sit inside one of these sets
    _within = {(_x, _y) for _m in _members for _n, _x in enumerate(_m) for _y in _m[_n + 1 :]}
    _explained = np.mean([
        (_x, _y) in _within
        for _x, _y in zip(guides.index.get_indexer(guides_together["guide_a"]), guides.index.get_indexer(guides_together["guide_b"]))
    ])
    _cross_replicate = int((_groups["replicate"].nunique() > 1).sum())
    _repeated_barcode = int((_groups["barcode"].nunique() < _groups.size()).sum())
    identical_sets = identical_sets.sort_values(["cells", "channels"], ascending=False, ignore_index=True)

    # gene UMIs of cells in one-channel groups, relative to the median cell with as many guides
    _umis = mdata["gene"].obs.loc[_obs.index, "total_counts"].to_numpy()
    _typical = pd.Series(_umis).groupby(_k).median()
    _at = _obs.index.get_indexer(_sets.index[_groups["channel_name"].transform("nunique").to_numpy() == 1])
    _depth = np.median(_umis[_at] / _typical.loc[_k[_at]].to_numpy())

    _observed = (identical_sets["channels"] == 1).mean()
    _chance = identical_sets["one_channel_by_chance"].mean()

    mo.vstack(
        [
            plot_one_channel(
                _observed,
                _chance,
                title=f"{_observed:.0%} of identical-guide groups sit in one channel",
                subtitle=f"Day 2: {len(identical_sets):,} groups, each of cells carrying exactly the same 3 or more guides",
            ),
            mo.md(
                f"**{len(identical_sets):,} groups** of cells carry exactly the same guide set: "
                f"{len(_sets) / len(_rows):.1%} of the {len(_rows):,} cells with 3 or more guides. Most groups are two "
                f"cells ({int((identical_sets['cells'] == 2).sum()):,}); the largest has {identical_sets['cells'].max()}. "
                f"{_explained:.0%} of the guide pairs found together in the previous step are two guides of one of "
                f"these sets. Removing all {len(_sets):,} of their cells would leave {len(_obs) - len(_sets):,} of the "
                f"{len(_obs):,} cells ({1 - len(_sets) / len(_obs):.1%}).\n\n"
                + (
                    "No group spans two replicates or reuses a barcode, so none is one library sequenced twice."
                    if _cross_replicate == 0 and _repeated_barcode == 0
                    else f"{_cross_replicate:,} groups span two replicates and {_repeated_barcode:,} reuse a barcode."
                )
                + "\n\n"
                f"- **Clones.** {1 - _observed:.0%} of groups have cells in two or more channels, which only clones can do.\n"
                f"- **One cell, several barcodes.** {_observed:.0%} sit in one channel, where channel sizes alone would "
                f"put {_chance:.0%}. Their cells carry {_depth:.2f}× the gene UMIs of cells with as many guides, "
                "consistent with each barcode holding part of one cell.\n\n"
                "For any single group, the data can't say which."
            ),
            mo.ui.table(identical_sets, selection=None),
        ]
    )
    return


@app.cell(hide_code=True)
def md_unassigned(mo):
    mo.md(r"""
    ### Guides in no cell

    The guides at zero, by type: those seen in some cell but never reaching the
    threshold, and those never seen at all.
    """)
    return


@app.cell
def unassigned_guides(assignment, cells_per_guide, guides, mdata, mo, np, pd):
    _counts = mdata["guide"].X
    _detected = np.bincount(_counts.indices, minlength=_counts.shape[1])
    _umis = np.asarray(_counts.sum(axis=0)).ravel()
    _never = cells_per_guide == 0

    _by_type = (
        pd.DataFrame(
            {
                "type": guides["type"].to_numpy(),
                "never_assigned": _never.to_numpy(),
                "never_seen": _detected == 0,
            }
        )
        .groupby("type", observed=True)
        .agg(guides=("never_assigned", "size"), never_assigned=("never_assigned", "sum"), never_seen=("never_seen", "sum"))
    )
    _unassigned = guides.loc[_never, ["type", "genomic_element", "intended_target_name"]].assign(
        cells_detected=_detected[_never.to_numpy()], total_umis=_umis[_never.to_numpy()]
    )

    mo.vstack(
        [
            mo.md(
                f"At **{assignment.threshold:.0f} UMIs**, **{int(_never.sum()):,}** of "
                f"{len(cells_per_guide):,} guides are assigned to no cell; "
                f"**{int((_detected == 0).sum()):,}** are never seen at all."
            ),
            _by_type,
            mo.ui.table(_unassigned.reset_index(names="guide_id")),
        ]
    )
    return


@app.cell(hide_code=True)
def md_cells_per_element(mo):
    mo.md(r"""
    ## Cells per targeted element

    An element is tested in every cell that carries any of its guides, so its
    cells are the union of its guides' cells: a cell with two guides of one element
    counts once. Positive controls and targeting elements together, each element
    keyed by its coordinates.
    """)
    return


@app.cell
def cells_per_element(assignment, guides, np, pd, plot_cells_per, sparse):
    _targeted = guides["intended_target_chr"].notna().to_numpy()
    _g = guides.loc[_targeted].assign(
        chr=lambda d: d["intended_target_chr"].astype(str),
        start=lambda d: d["intended_target_start"].astype("int64"),
        end=lambda d: d["intended_target_end"].astype("int64"),
    )
    _codes, _ids = pd.factorize(_g["chr"] + ":" + _g["start"].astype(str) + "-" + _g["end"].astype(str))
    # elements × guides: 1 where the guide targets the element
    _member = sparse.csr_matrix(
        (np.ones(len(_codes), dtype=np.int32), (_codes, np.flatnonzero(_targeted))),
        shape=(len(_ids), len(guides)),
    )
    # elements × cells: a cell carrying two guides of one element counts once
    element_cells = (_member @ assignment.assigned.astype(np.int32)).tocsr()
    element_cells.sort_indices()

    elements = (
        _g[["chr", "start", "end", "type", "genomic_element", "element_location", "intended_target_name", "putative_target_genes"]]
        .groupby(_codes, observed=True)
        .first()
        .set_axis(pd.Index(_ids, name="element"))
        .assign(guides=np.diff(_member.indptr), cells=np.diff(element_cells.indptr))
    )
    plot_cells_per(
        elements["cells"],
        unit="element",
        bin_width=20,
        cap=1000,
        title=f"The median element is tested in {elements['cells'].median():,.0f} cells",
        subtitle=(
            f"Day 2: cells carrying any guide of each of the {len(elements):,} targeted elements "
            f"({assignment.threshold:.0f}-UMI threshold)"
        ),
    )
    return element_cells, elements


@app.cell(hide_code=True)
def md_elements_by_guides(mo):
    mo.md(r"""
    ### Elements by guides per element

    Where the elements with few cells come from.
    """)
    return


@app.cell
def elements_by_guides(elements, np, pd):
    _bins = pd.cut(elements["guides"], [0, 1, 5, 10, 15, np.inf], labels=["1", "2–5", "6–10", "11–15", "16+"])
    elements.groupby(_bins, observed=True).agg(
        elements=("cells", "size"),
        cells_median=("cells", "median"),
        under_100_cells=("cells", lambda c: int((c < 100).sum())),
    ).rename_axis("guides per element")
    return


@app.cell(hide_code=True)
def md_positive_controls(mo):
    mo.md(r"""
    ## Pairs: positive controls

    pysceptre tests **pairs**: an element (`grna_target`, the cells carrying any of
    its guides) against a gene (`response_id`). Positive controls are the pairs
    whose effect we already expect: each positive-control element of the design,
    its promoters and validated enhancers, against the gene the design names for
    it. Only readout genes are paired: TAP-seq measures nothing else.
    """)
    return


@app.cell
def positive_controls(elements, genes, mo, pair_table):
    readout_genes = (
        genes.loc[genes["readout"], ["symbol", "gene_chr", "gene_start", "gene_end", "gene_tss", "gene_strand"]]
        .astype({"gene_chr": str, "gene_strand": str, "gene_tss": "int64"})
        .assign(gene_id=lambda d: d.index.str.replace(r"\.\d+$", "", regex=True))
        .rename_axis("response_id")
    )

    # each positive-control element with each gene the design names for it
    _named = (
        elements.loc[elements["type"] == "positive control", "putative_target_genes"]
        .astype(str)
        .str.split(",")
        .explode()
        .rename("gene_id")
        .rename_axis("grna_target")
        .reset_index()
    )
    _measured = _named.merge(readout_genes["gene_id"].reset_index(), on="gene_id")
    positive_controls = pair_table(_measured, "positive_control", elements, readout_genes)

    mo.vstack(
        [
            mo.md(
                f"The design names **{len(_named)} pairs** for its {_named['grna_target'].nunique()} "
                f"positive-control elements. **{len(positive_controls)}** have a readout gene and are the "
                f"positive controls; the other {len(_named) - len(positive_controls)} name a gene TAP-seq "
                "does not measure."
            ),
            positive_controls.groupby("genomic_element", observed=True).agg(
                pairs=("grna_target", "size"),
                elements=("grna_target", "nunique"),
                genes=("response_id", "nunique"),
            ),
            mo.ui.table(positive_controls),
        ]
    )
    return positive_controls, readout_genes


@app.cell(hide_code=True)
def md_host_genes(mo):
    mo.md(r"""
    ## Pairs: host genes

    An element inside a gene's body is paired with that gene. CRISPRi there is
    expected to have a negative effect on the gene, so these pairs are positive controls: they go to the
    power check, not to discovery. They keep their own type, `host_gene`, so they
    are reported apart from the pairs the design names. A pair the design already
    names stays a `positive_control`.
    """)
    return


@app.cell
def host_gene_pairs(
    elements,
    mo,
    pair_table,
    pd,
    positive_controls,
    readout_genes,
):
    # every element with every readout gene on its chromosome, kept where it overlaps the gene's body
    _by_chr = (
        elements.rename_axis("grna_target")
        .reset_index()
        .merge(readout_genes.reset_index(), left_on="chr", right_on="gene_chr")
    )
    _inside = _by_chr[(_by_chr["start"] < _by_chr["gene_end"]) & (_by_chr["end"] > _by_chr["gene_start"])]

    # a pair the design already names stays a positive control
    _designed = pd.MultiIndex.from_frame(_inside[["grna_target", "response_id"]]).isin(
        pd.MultiIndex.from_frame(positive_controls[["grna_target", "response_id"]])
    )
    host_gene_pairs = pair_table(_inside[~_designed], "host_gene", elements, readout_genes)

    _genic = (elements["type"] == "targeting") & (elements["element_location"] == "genic")
    _hosted = elements.index.isin(_inside["grna_target"])

    mo.vstack(
        [
            mo.md(
                f"**{len(_inside):,}** pairs of an element and the readout gene it sits in; "
                f"{int(_designed.sum())} are already positive controls, which leaves "
                f"**{len(host_gene_pairs):,} host-gene pairs**: "
                f"{host_gene_pairs['grna_target'].nunique():,} elements in "
                f"{host_gene_pairs['response_id'].nunique()} readout genes. "
                f"{'All' if elements.loc[_hosted, 'element_location'].eq('genic').all() else 'Not all'} "
                "of these elements are labelled genic in the guide table; "
                f"{int((_genic & ~_hosted).sum()):,} genic targeting elements sit only in genes "
                "TAP-seq does not measure, so they have no host-gene pair."
            ),
            host_gene_pairs.groupby("genomic_element", observed=True).agg(
                pairs=("grna_target", "size"),
                elements=("grna_target", "nunique"),
                genes=("response_id", "nunique"),
            ),
            mo.ui.table(host_gene_pairs),
        ]
    )
    return (host_gene_pairs,)


@app.cell(hide_code=True)
def md_discovery(mo):
    mo.md(r"""
    ## Pairs: discovery

    Every element is tested against each readout gene whose TSS lies within 1 Mb of
    the element's midpoint, unless the pair is already a positive control or a
    host-gene pair. As in the DC-TAP-seq paper, every discovery pair enters the
    multiple-testing correction except a **self promoter**: a promoter element
    against its own gene, one of the genes `putative_target_genes` names for it.
    Self promoters are tested, but `in_fdr` is `False` for them. A promoter element
    against another gene stays in.
    """)
    return


@app.cell
def discovery_pairs(
    elements,
    host_gene_pairs,
    mo,
    np,
    pair_table,
    pd,
    positive_controls,
    readout_genes,
):
    MAX_DISTANCE = 1_000_000  # bases from the element's midpoint to the gene's TSS

    # every element with every readout gene on its chromosome, kept where the TSS is within MAX_DISTANCE
    _by_chr = (
        elements.rename_axis("grna_target")
        .reset_index()
        .merge(readout_genes.reset_index(), left_on="chr", right_on="gene_chr")
    )
    _near = _by_chr[((_by_chr["start"] + _by_chr["end"]) // 2 - _by_chr["gene_tss"]).abs() <= MAX_DISTANCE]

    # a pair already tested as a positive control or a host gene is not tested again
    _taken = pd.MultiIndex.from_frame(pd.concat([positive_controls, host_gene_pairs])[["grna_target", "response_id"]])
    _new = ~pd.MultiIndex.from_frame(_near[["grna_target", "response_id"]]).isin(_taken)
    discovery_pairs = pair_table(_near[_new], "discovery", elements, readout_genes)

    # as in the DC-TAP-seq paper, every pair enters the multiple-testing correction but a self promoter:
    # a promoter element against its own gene, one of the genes putative_target_genes names for it
    _own_genes = (
        discovery_pairs["grna_target"].map(elements["putative_target_genes"].astype("string")).fillna("").str.split(",")
    )
    _gene = discovery_pairs["response_id"].str.replace(r"\.\d+$", "", regex=True)
    _self_promoter = (discovery_pairs["genomic_element"] == "promoter").to_numpy() & np.array(
        [_g in _own for _g, _own in zip(_gene, _own_genes)]
    )
    discovery_pairs["in_fdr"] = ~_self_promoter

    _paired = pd.concat([positive_controls, host_gene_pairs, discovery_pairs])["grna_target"]
    _unpaired = int((~elements.index.isin(_paired)).sum())

    mo.vstack(
        [
            mo.md(
                f"**{len(discovery_pairs):,} discovery pairs**: "
                f"{discovery_pairs['grna_target'].nunique():,} elements against "
                f"{discovery_pairs['response_id'].nunique()} readout genes. "
                f"**{int(discovery_pairs['in_fdr'].sum()):,}** enter the multiple-testing correction; "
                f"{int(_self_promoter.sum()):,} self promoters are tested but do not. "
                f"{int(_near.shape[0] - len(discovery_pairs)):,} pairs within 1 Mb are already positive "
                f"controls or host-gene pairs, and {_unpaired:,} elements have no pair of any type."
            ),
            discovery_pairs.groupby("genomic_element", observed=True).agg(
                pairs=("grna_target", "size"),
                elements=("grna_target", "nunique"),
                genes=("response_id", "nunique"),
                in_fdr=("in_fdr", "sum"),
            ),
            mo.ui.table(discovery_pairs),
        ]
    )
    return (discovery_pairs,)


@app.cell(hide_code=True)
def md_pysceptre_inputs(mo):
    mo.md(r"""
    ## Inputs for pysceptre

    pysceptre takes plain arrays, not the MuData object:

    - `response_matrix`: genes × cells counts, here the readout genes only;
    - `gene_ids` and `cell_ids`: its row and column names, the genes and the
      barcodes. pysceptre is given `gene_ids` but never the barcodes: it works on
      positions, and `cell_ids` maps them back;
    - `grna_target_cells`: for each element, the positions of the cells carrying any
      of its guides;
    - `pairs`: the element–gene pairs from the three steps above.

    Cells are referred to by position, so the gene and guide matrices must list
    them in the same order.
    """)
    return


@app.cell
def pysceptre_inputs(
    discovery_pairs,
    element_cells,
    elements,
    host_gene_pairs,
    mdata,
    mo,
    np,
    pd,
    positive_controls,
    readout_genes,
):
    # pysceptre refers to cells by position, so both modalities must list them in the same order
    assert mdata["gene"].obs_names.equals(mdata["guide"].obs_names)

    cell_ids = mdata["gene"].obs_names  # barcodes: column i of response_matrix is cell_ids[i]
    gene_ids = list(readout_genes.index)  # genes: row i of response_matrix is gene_ids[i]
    response_matrix = mdata["gene"][:, gene_ids].X.T.tocsr()  # readout genes × cells

    grna_target_cells = {
        _target: element_cells.indices[element_cells.indptr[_i] : element_cells.indptr[_i + 1]]
        for _i, _target in enumerate(elements.index)
    }

    pairs = pd.concat([positive_controls, host_gene_pairs, discovery_pairs], ignore_index=True)
    pairs["in_fdr"] = pairs["in_fdr"].astype("boolean")  # set for discovery pairs only
    assert not pairs.duplicated(["grna_target", "response_id"]).any()

    mo.vstack(
        [
            mo.md(
                f"- `response_matrix`: {response_matrix.shape[0]} readout genes × "
                f"{response_matrix.shape[1]:,} cells, {response_matrix.nnz:,} nonzero\n"
                f"- `gene_ids`: its {len(gene_ids)} rows, versioned Ensembl ids (symbols in `readout_genes`)\n"
                f"- `cell_ids`: its {len(cell_ids):,} columns, the barcodes, e.g. `{cell_ids[0]}`\n"
                f"- `grna_target_cells`: {len(grna_target_cells):,} elements, a median of "
                f"{np.median([len(_c) for _c in grna_target_cells.values()]):,.0f} cells each, as "
                "positions in `cell_ids`\n"
                f"- `pairs`: {len(pairs):,}, no element–gene pair twice; `in_fdr` is set for the "
                f"discovery pairs only, {int(pairs['in_fdr'].sum()):,} of which enter the "
                "multiple-testing correction"
            ),
            pairs.groupby("pair_type", sort=False).agg(
                pairs=("grna_target", "size"),
                elements=("grna_target", "nunique"),
                genes=("response_id", "nunique"),
            ),
        ]
    )
    return cell_ids, gene_ids, grna_target_cells, pairs, response_matrix


@app.cell(hide_code=True)
def md_pairwise_qc(mo):
    mo.md(r"""
    ## Pairwise QC

    A pair can be tested only if its gene is detected in enough cells on both
    sides: at least 20 **treated** cells (carrying any of the element's guides) and
    90,000 **control** cells (all the other cells) with a nonzero count. sceptre's
    defaults are 7 and 7. These are day 0's cuts: there the readout genes split at
    about 100,000 nonzero control cells, and the control cut sat in the gap. Day 2
    has no such gap, as the scatter below shows; we keep the same cuts on every
    day, for consistency. pysceptre's power check applies this itself; discovery expects pairs that
    already pass, so we apply it here, with the same function, to every pair.
    """)
    return


@app.cell
def pairwise_qc(
    gene_ids,
    grna_target_cells,
    mo,
    np,
    pairs,
    pd,
    pysceptre,
    response_matrix,
):
    N_NONZERO_TRT, N_NONZERO_CNTRL = 20, 90_000

    # the function pysceptre's power check uses; control cells are all the others (high MOI)
    pairs_qc = pysceptre.pipeline.power.annotate_pairwise_qc(
        pairs,
        response_matrix,
        gene_ids,
        grna_target_cells,
        response_matrix.shape[1],
        n_nonzero_trt_thresh=N_NONZERO_TRT,
        n_nonzero_cntrl_thresh=N_NONZERO_CNTRL,
    )

    _summary = pairs_qc.groupby("pair_type", sort=False).agg(
        pairs=("pass_qc", "size"),
        pass_qc=("pass_qc", "sum"),
    )
    _summary["fail_qc"] = _summary["pairs"] - _summary["pass_qc"]

    # where the failures come from: the share of cells in which each gene is detected
    _detected = pd.Series(
        100 * np.diff(response_matrix.indptr) / response_matrix.shape[1], index=gene_ids, name="cells_detected_pct"
    )
    _by_gene = (
        pairs_qc.groupby(["response_id", "gene_symbol"])
        .agg(pairs=("pass_qc", "size"), fail_qc=("pass_qc", lambda _p: int((~_p).sum())))
        .reset_index("gene_symbol")
        .join(_detected)
    )
    _failing = _by_gene["fail_qc"] > 0
    _fails = pairs_qc.loc[~pairs_qc["pass_qc"]]

    mo.vstack(
        [
            mo.md(
                f"**{int(pairs_qc['pass_qc'].sum()):,} of {len(pairs_qc):,} pairs pass**; "
                f"{len(_fails):,} fail. "
                f"{int((pairs_qc['pass_qc'] & pairs_qc['in_fdr'].fillna(False)).sum()):,} discovery pairs "
                "that pass enter the multiple-testing correction."
            ),
            _summary,
            mo.md(
                f"{'Every' if (_fails['n_nonzero_trt'] < N_NONZERO_TRT).all() else 'Not every'} failure is "
                f"on the treated side, and they fall on {int(_failing.sum())} of {len(gene_ids)} genes. "
                f"Those genes are detected in a median "
                f"{_by_gene.loc[_failing, 'cells_detected_pct'].median():.1f}% of cells, against "
                f"{_by_gene.loc[~_failing, 'cells_detected_pct'].median():.0f}% for the genes whose pairs "
                "all pass."
            ),
            mo.ui.table(_by_gene.loc[_failing].sort_values("fail_qc", ascending=False).round(2)),
        ]
    )
    return N_NONZERO_CNTRL, N_NONZERO_TRT, pairs_qc


@app.cell(hide_code=True)
def md_pairwise_qc_scatter(mo):
    mo.md(r"""
    ### Treated against control cells, pair by pair

    Each pair's two counts that the QC cuts on: the cells carrying the element's
    guides (treated) and all the other cells (control) that have a nonzero count of
    the pair's gene. A pair needs enough of each.
    """)
    return


@app.cell
def pairwise_qc_scatter(
    N_NONZERO_CNTRL,
    N_NONZERO_TRT,
    mo,
    pairs_qc,
    plot_pairwise_qc,
):
    _fail = int((~pairs_qc["pass_qc"]).sum())
    _on_control = int((pairs_qc["n_nonzero_cntrl"] < N_NONZERO_CNTRL).sum())
    mo.vstack(
        [
            plot_pairwise_qc(
                pairs_qc,
                N_NONZERO_TRT,
                N_NONZERO_CNTRL,
                title=f"{_fail:,} pairs fail, "
                + (f"{_on_control:,} of them on the control cut" if _on_control else "all on the treated cut"),
                subtitle="Day 2: cells with a nonzero count of the pair's gene; log scales",
            ),
            mo.md(
                "One dot per pair. Treated: cells carrying any guide of the pair's element. Control: all the other cells. "
                "Each axis counts only the cells with a nonzero count of the pair's gene. "
                f"Dashed lines: the cuts, {N_NONZERO_TRT} treated and {N_NONZERO_CNTRL} control cells; "
                "a pair below either line fails."
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def md_pairs_by_expression(mo):
    mo.md(r"""
    ### Pairs by the gene's expression

    How many pairs fall at each expression level of their gene, and how many of them
    pass. Expression is the gene's mean UMIs per cell over every cell kept, in four
    bins.
    """)
    return


@app.cell
def pairs_by_expression(gene_ids, mo, np, pairs_qc, pd, response_matrix):
    # each readout gene's mean UMIs per cell, over every cell kept
    _mean = pd.Series(np.asarray(response_matrix.sum(axis=1)).ravel() / response_matrix.shape[1], index=gene_ids)
    _bins = pd.cut(_mean, [0, 0.1, 1, 10, np.inf], right=False, labels=["< 0.1", "0.1–1", "1–10", "≥ 10"])

    _pairs = pairs_qc.assign(
        expression=pairs_qc["response_id"].map(_bins).astype(_bins.dtype),
        pair_type=pd.Categorical(pairs_qc["pair_type"], categories=pairs_qc["pair_type"].unique()),
    )
    pairs_by_expression = _pairs.groupby(["expression", "pair_type"], observed=False).agg(
        genes=("response_id", "nunique"), pairs=("pass_qc", "size"), pass_qc=("pass_qc", "sum")
    )
    pairs_by_expression["fail_qc"] = pairs_by_expression["pairs"] - pairs_by_expression["pass_qc"]

    mo.vstack(
        [
            mo.md(
                f"Of the {len(gene_ids)} readout genes, by mean UMIs per cell: "
                + ", ".join(f"**{_n}** at {_b}" for _b, _n in _bins.value_counts(sort=False).items())
                + "."
            ),
            pairs_by_expression,
        ]
    )
    return


@app.cell
def md_covariates(mo):
    mo.md(r"""
    ## Covariates

    pysceptre tests each pair with a regression, so whatever else moves a gene's
    count in a cell goes in as a covariate. We use sceptre's own high-MOI design,
    each count as a log:

    - **sequencing depth**: the cell's gene UMIs and genes detected, over every
      gene, not only the readout;
    - **guide load**: the cell's guide UMIs and guides detected (every guide with a
      UMI, assigned or not).

    To these we add the **channel**. Each channel is its own library and sequencing
    run, and channels sit inside replicates, so the channel term also covers the
    differences between replicates. sceptre's automatic formula drops any factor
    with 15 or more levels, so we add the channel by hand: one column per channel
    but the first, which the intercept stands for.

    pysceptre takes the design as a plain matrix, `covariate_matrix`: one row per
    cell, in `cell_ids` order, intercept included.
    """)
    return


@app.cell
def covariates(cell_ids, mdata, mo, np, pd):
    _genes = mdata["gene"].obs.loc[cell_ids]
    _guides = mdata["guide"].obs.loc[cell_ids]

    # sceptre's high-MOI covariates, in its order; every cell has at least one of each, so no log(x + 1)
    covariates = pd.DataFrame(
        {
            "log_genes_detected": np.log(_genes["n_genes_by_counts"].to_numpy(dtype=float)),
            "log_gene_umis": np.log(_genes["total_counts"].to_numpy(dtype=float)),
            "log_guides_detected": np.log(_guides["n_guides_detected"].to_numpy(dtype=float)),
            "log_guide_umis": np.log(_guides["total_guide_umis"].to_numpy(dtype=float)),
            "channel": mdata.obs.loc[cell_ids, "channel_name"].astype(str).to_numpy(),
        },
        index=cell_ids,
    )
    # one column per channel but the first, which the intercept stands for
    _channels = pd.get_dummies(covariates["channel"], prefix="channel", drop_first=True, dtype=float)
    _logs = covariates.drop(columns="channel")

    covariate_matrix = np.column_stack([np.ones(len(cell_ids)), _logs.to_numpy(), _channels.to_numpy()])
    covariate_names = ["intercept", *_logs.columns, *_channels.columns]
    _rank = np.linalg.matrix_rank(covariate_matrix)

    mo.vstack(
        [
            mo.md(
                f"`covariate_matrix`: {covariate_matrix.shape[0]:,} cells × {covariate_matrix.shape[1]} columns, "
                f"the intercept, {_logs.shape[1]} log counts and {_channels.shape[1]} channel columns "
                f"(reference channel `{sorted(covariates['channel'].unique())[0]}`). "
                + (
                    "Full rank: no column is a combination of the others."
                    if _rank == covariate_matrix.shape[1]
                    else f"**Rank {_rank} of {covariate_matrix.shape[1]}**: some column is a combination of others."
                )
                + "\n\nMedian per cell, by replicate:"
            ),
            np.exp(_logs)
            .set_axis(["genes detected", "gene UMIs", "guides detected", "guide UMIs"], axis=1)
            .groupby(mdata.obs.loc[cell_ids, "replicate"].astype(str).to_numpy())
            .median()
            .round()
            .astype(int)
            .rename_axis("replicate"),
        ]
    )
    return


@app.cell(hide_code=True)
def md_save(mo):
    mo.md(r"""
    ## Save the cells kept

    The cells kept, as loaded: `gene` (UMIs per gene, with `gene_tss` and
    `gene_strand`) and `guide` (UMIs per guide and the `guide_assignment` layer),
    with `pass_qc`, `pass_limits` and `keep` in `obs`. `uns["cell_limits"]` records
    the three limits, `uns["source"]` the day object. Writing takes minutes, so it
    runs only when the button is pressed.
    """)
    return


@app.cell
def save_button(mo):
    save_button = mo.ui.run_button(label="Write the h5mu")
    save_button
    return (save_button,)


@app.cell
def save_h5mu(
    MIN_GENE_UMIS,
    MIN_READOUT_PCT,
    MIN_RIBO_PCT,
    data_path,
    mdata,
    mo,
    save_button,
):
    mo.stop(not save_button.value, mo.md("Press the button to write the file."))

    OUT = mo.notebook_dir().parents[1] / "results" / "2026-11-23_day_objects_qc" / f"{data_path.stem}_pass_qc.h5mu"
    OUT.parent.mkdir(parents=True, exist_ok=True)
    mdata.uns["cell_limits"] = {
        "min_gene_umis": MIN_GENE_UMIS,
        "min_readout_pct": MIN_READOUT_PCT,
        "min_ribo_pct": MIN_RIBO_PCT,
    }
    mdata.uns["source"] = str(data_path)
    _partial = OUT.with_suffix(".partial.h5mu")
    mdata.write_h5mu(_partial, compression="gzip")
    _partial.replace(OUT)  # only a finished write takes the final name
    mo.md(f"Wrote **{mdata.n_obs:,}** cells to `{OUT}` ({OUT.stat().st_size / 2**30:.1f} GiB).")
    return


@app.cell(hide_code=True)
def md_helpers(mo):
    mo.md(r"""
    ## Helpers

    The functions the steps above call. marimo runs cells in the order their
    variables require, not in page order, so they can live down here.
    """)
    return


@app.cell
def knee_plot_fn(matplotlib, np, plt):
    def plot_knees(curves, cut=None, ylabel="UMIs per barcode", title="", subtitle=""):
        """Each channel's barcodes ranked by UMIs (left) and their geometric mean ± 1 SD (right).

        The mean is taken over log10 UMIs, only at ranks every channel reaches.
        """
        blue, gridline = "#006BA2", "#B7C6CF"
        fig, (left, right) = plt.subplots(
            1, 2, figsize=(7.5, 3.8), sharey=True, gridspec_kw={"wspace": 0.14}
        )
        fig.subplots_adjust(left=0.11, right=0.95, bottom=0.14, top=0.74)

        grid = np.unique(np.geomspace(1, curves["rank"].max(), 200).round().astype(int))
        logs = []
        for _, channel in curves.groupby("measurement_set", observed=True):
            rank = channel["rank"].to_numpy(float)
            umis = channel["umis"].to_numpy(float)
            left.plot(rank, umis, color="0.55", linewidth=0.6, alpha=0.5)
            reach = grid <= rank.max()
            values = np.full(len(grid), np.nan)
            values[reach] = np.interp(np.log10(grid[reach]), np.log10(rank), np.log10(umis))
            logs.append(values)
        logs = np.vstack(logs)
        shared = ~np.isnan(logs).any(axis=0)
        mean = logs[:, shared].mean(axis=0)
        sd = logs[:, shared].std(axis=0, ddof=1)
        right.fill_between(
            grid[shared], 10 ** (mean - sd), 10 ** (mean + sd), color=blue, alpha=0.2, linewidth=0
        )
        right.plot(grid[shared], 10**mean, color=blue, linewidth=1.5)

        plain = matplotlib.ticker.FuncFormatter(lambda value, _: f"{value:,.0f}")
        for ax in (left, right):
            ax.set(xscale="log", yscale="log", xlabel="Barcode rank")
            ax.xaxis.set_major_locator(matplotlib.ticker.LogLocator(base=100))
            ax.xaxis.set_major_formatter(plain)
            ax.yaxis.set_major_formatter(plain)
            ax.minorticks_off()
            ax.grid(axis="y", linewidth=0.5, color=gridline)
            ax.grid(axis="x", visible=False)
            ax.spines[["top", "right", "left"]].set_visible(False)
            ax.tick_params(axis="y", length=0)
            if cut:
                ax.axhline(cut, color="0.1", linewidth=1, linestyle=(0, (3, 1)))
        left.set_ylabel(ylabel)
        left.set_title(f"{len(logs)} channels", loc="left", fontsize=9, color="0.3")
        right.set_title("Geometric mean ± 1 SD", loc="left", fontsize=9, color="0.3")
        if cut:
            right.annotate(
                f"cut: {cut:,.0f} UMIs",
                (grid[0], cut),
                xytext=(2, 3),
                textcoords="offset points",
                fontsize=8,
                color="0.1",
            )
        fig.text(0.11, 0.94, title, fontsize=13, fontweight="bold", ha="left")
        fig.text(0.11, 0.885, subtitle, fontsize=9.5, color="0.3", ha="left")
        return fig

    return (plot_knees,)


@app.cell
def share_plot_fn(matplotlib, np, plt):
    def plot_share(share, mark=None, xlabel="", title="", subtitle=""):
        """Histogram of a per-cell share (0 to 1), with its median and an optional reference line."""
        blue, gridline = "#006BA2", "#B7C6CF"
        fig, ax = plt.subplots(figsize=(7.5, 3.6))
        fig.subplots_adjust(left=0.11, right=0.95, bottom=0.15, top=0.76)
        percent = 100 * np.asarray(share)
        heights, _, _ = ax.hist(percent, bins=np.arange(0, 101, 1), color=blue, linewidth=0)
        ax.set(xlim=(0, 100), ylim=(0, heights.max() * 1.2), xlabel=xlabel, ylabel="Cells")

        def label(x, text, side):
            ax.axvline(x, ymax=1 / 1.2 + 0.08, color="0.1", linewidth=1,
                       linestyle="-" if side > 0 else (0, (3, 1)))
            ax.annotate(text, (x, 1), xycoords=("data", "axes fraction"), xytext=(4 * side, -9),
                        textcoords="offset points", fontsize=8, color="0.1",
                        ha="left" if side > 0 else "right")

        median = np.median(percent)
        label(median, f"median {median:.0f}%", +1)
        if mark is not None:
            below = int((percent < 100 * mark).sum())
            label(100 * mark, f"{below:,} cells below {100 * mark:.0f}%", -1)

        ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
        ax.grid(axis="y", linewidth=0.5, color=gridline)
        ax.grid(axis="x", visible=False)
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.tick_params(axis="y", length=0)
        fig.text(0.11, 0.93, title, fontsize=13, fontweight="bold", ha="left")
        fig.text(0.11, 0.875, subtitle, fontsize=9.5, color="0.3", ha="left")
        return fig

    return (plot_share,)


@app.cell
def share_by_plot_fn(np, pd, plt):
    def plot_share_by(share, groups, label="Batch", row_height=0.85, xlim=(0, 100), xlabel="", title="", subtitle=""):
        """One histogram row per group of a per-cell share (0 to 1), each as % of its own cells."""
        blue, gridline = "#006BA2", "#B7C6CF"
        groups = pd.Series(np.asarray(groups), index=share.index)
        names = sorted(groups.unique(), key=str)
        gap = 0.3  # inches between rows, room for each row's title
        height = 1.5 + (row_height + gap) * len(names)
        bottom, top = 0.55 / height, 1 - 0.95 / height
        axes_height = (height * (top - bottom) - gap * (len(names) - 1)) / len(names)
        fig, axes = plt.subplots(len(names), 1, figsize=(7.5, height), sharex=True, sharey=True)
        fig.subplots_adjust(left=0.11, right=0.95, bottom=bottom, top=top, hspace=gap / axes_height)
        bins = np.arange(xlim[0], xlim[1] + 1, 1)
        tallest = 0
        for ax, name in zip(axes, names):
            percent = 100 * share[groups == name].to_numpy()
            heights, _, _ = ax.hist(
                percent, bins=bins, weights=np.full(len(percent), 100 / len(percent)),
                color=blue, linewidth=0,
            )
            tallest = max(tallest, heights.max())
            median = np.median(percent)
            ax.axvline(median, color="0.1", linewidth=1)
            ax.annotate(
                f"median {median:.0f}%", (median, 1), xycoords=("data", "axes fraction"),
                xytext=(4, -9), textcoords="offset points", fontsize=8, color="0.1",
            )
            ax.set_title(f"{label} {name} · {len(percent):,} cells".strip(), loc="left", fontsize=9, color="0.3")
            ax.grid(axis="y", linewidth=0.5, color=gridline)
            ax.grid(axis="x", visible=False)
            ax.spines[["top", "right", "left"]].set_visible(False)
            ax.tick_params(axis="y", length=0)
        axes[0].set(xlim=xlim, ylim=(0, tallest * 1.3))
        axes[-1].set_xlabel(xlabel)
        fig.text(0.02, 0.5, "Cells, % of the row's cells", rotation=90, va="center", fontsize=10)
        fig.text(0.11, 1 - 0.3 / height, title, fontsize=13, fontweight="bold", ha="left")
        fig.text(0.11, 1 - 0.52 / height, subtitle, fontsize=9.5, color="0.3", ha="left")
        return fig

    return (plot_share_by,)


@app.cell
def r_qc_plot_fns(headline, house_style, jscatter, matplotlib, mo, np, plt):
    def plot_hist_median(values, bins=50, xlabel="", title="", subtitle=""):
        """Histogram of a per-cell value, in cells, with its median."""
        values = np.asarray(values, dtype=float)
        fig, ax = plt.subplots(figsize=(7.5, 3.6))
        fig.subplots_adjust(left=0.11, right=0.95, bottom=0.15, top=0.76)
        ax.hist(values, bins=bins, color="#006BA2", linewidth=0)
        median = np.median(values)
        ax.axvline(median, color="0.1", linewidth=1)
        ax.annotate(f"median {median:.2f}%", (median, 1), xycoords=("data", "axes fraction"), xytext=(4, -9),
                    textcoords="offset points", fontsize=8, color="0.1")
        ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
        ax.set(xlabel=xlabel, ylabel="Cells")
        ax.tick_params(labelsize=8)
        house_style(ax)
        headline(fig, title, subtitle)
        return fig


    def plot_upm_box(upm, readout, title="", subtitle=""):
        """Mean UPM plus 1 of the readout genes and of all other genes, as box plots on a log scale."""
        readout = np.asarray(readout, dtype=bool)
        upm = np.asarray(upm, dtype=float)
        groups = {"Readout genes": upm[readout], "Other genes": upm[~readout]}
        fig, ax = plt.subplots(figsize=(7.5, 3.8))
        fig.subplots_adjust(left=0.11, right=0.95, bottom=0.15, top=0.76)
        ax.boxplot(
            [values + 1 for values in groups.values()],
            tick_labels=[f"{name}\n{len(values):,} genes" for name, values in groups.items()],
            widths=0.5,
            patch_artist=True,
            boxprops=dict(facecolor="#98DAFF", edgecolor="#006BA2"),
            whiskerprops=dict(color="#006BA2"),
            capprops=dict(color="#006BA2"),
            medianprops=dict(color="#A81829", linewidth=2),
            flierprops=dict(marker=".", markersize=2, markerfacecolor="#006BA2", markeredgewidth=0, alpha=0.3),
        )
        ax.set_yscale("log")
        ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:,g}"))
        ax.set(ylabel="Mean UPM + 1")
        ax.tick_params(axis="x", length=0)
        ax.tick_params(labelsize=8)
        house_style(ax)
        headline(fig, title, subtitle)
        return fig


    def plot_marker_strip(upm, classes, order, title="", subtitle="", seed=0):
        """Each gene's mean UPM by class, spread sideways, on a log scale, with the class median."""
        rng = np.random.default_rng(seed)
        upm, classes = np.asarray(upm, dtype=float), np.asarray(classes)
        fig, ax = plt.subplots(figsize=(7.5, 3.8))
        fig.subplots_adjust(left=0.11, right=0.95, bottom=0.15, top=0.76)
        for i, name in enumerate(order):
            values = upm[classes == name]
            shown = values[values > 0]
            ax.scatter(i + rng.uniform(-0.2, 0.2, len(shown)), shown, s=16, color="#006BA2", alpha=0.8, linewidths=0)
            median = np.median(values)
            if median > 0:
                ax.hlines(median, i - 0.3, i + 0.3, color="#A81829", linewidth=2)
                ax.annotate(f"{median:,.1f}", (i + 0.32, median), va="center", fontsize=8, color="#A81829")
        ax.set_xticks(range(len(order)), [f"{name}\n{int((classes == name).sum())} genes" for name in order])
        ax.set_yscale("log")
        ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:,g}"))
        ax.set(xlim=(-0.6, len(order) - 0.4), ylabel="Mean UPM")
        ax.tick_params(axis="x", length=0)
        ax.tick_params(labelsize=8)
        house_style(ax)
        headline(fig, title, subtitle)
        return fig


    def threshold_scatter(data, x, y, x_cut=None, y_cut=None, xlabel="", ylabel="", scale="linear", tooltip=("x", "y"),
                          line=None, title="", subtitle="", caption=""):
        """An interactive jscatter plot of two columns of `data` under a headline: lines at the cuts given, and a red line
        through the (x, y) vertices of `line`."""
        scatter = jscatter.Scatter(x=x, y=y, data=data, x_scale=scale, y_scale=scale, height=420)
        scatter.color(default="#006BA2")
        scatter.axes(True, grid=True, labels=[xlabel, ylabel])
        annotations = [jscatter.VLine(x_cut, line_color="#1A1A1A")] if x_cut is not None else []
        if y_cut is not None:
            annotations.append(jscatter.HLine(y_cut, line_color="#1A1A1A"))
        if line is not None:
            annotations.append(jscatter.Line(list(line), line_color="#A81829", line_width=2))
        if annotations:
            scatter.annotations(annotations)
        scatter.tooltip(True, properties=list(tooltip))
        return mo.vstack(
            [
                mo.md(f"**{title}**<br><span style='color:#595959'>{subtitle}</span>"),
                mo.ui.anywidget(scatter.widget),
                mo.md(caption),
            ]
        )

    return plot_hist_median, plot_marker_strip, plot_upm_box, threshold_scatter


@app.cell
def counts_plot_fn(np, plt):
    def plot_counts(counts, cap=20, xlabel="", title="", subtitle=""):
        """Bar per whole-number value, as % of cells; values above `cap` share one last bar."""
        blue, gridline = "#006BA2", "#B7C6CF"
        counts = np.asarray(counts)
        values = np.arange(counts.min(), cap + 2)
        tally = np.bincount(np.minimum(counts, cap + 1), minlength=cap + 2)[values]
        percent = 100 * tally / len(counts)
        fig, ax = plt.subplots(figsize=(7.5, 3.6))
        fig.subplots_adjust(left=0.11, right=0.95, bottom=0.15, top=0.76)
        ax.bar(values, percent, width=0.8, color=blue, linewidth=0)
        ax.set_xticks(values, [str(v) for v in values[:-1]] + [f"{cap + 1}+"])
        ax.set(xlim=(values[0] - 0.6, values[-1] + 0.6), xlabel=xlabel, ylabel="Cells, %")
        ax.grid(axis="y", linewidth=0.5, color=gridline)
        ax.grid(axis="x", visible=False)
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.tick_params(axis="y", length=0)
        ax.tick_params(axis="x", length=0, labelsize=8)
        fig.text(0.11, 0.93, title, fontsize=13, fontweight="bold", ha="left")
        fig.text(0.11, 0.875, subtitle, fontsize=9.5, color="0.3", ha="left")
        return fig

    return (plot_counts,)


@app.cell
def sweep_fns(np, plt):
    def guides_at(counts, threshold):
        """Guides per cell with at least `threshold` UMIs, from a cells × guides CSR matrix."""
        above = np.concatenate([[0], np.cumsum(counts.data >= threshold)])
        return np.diff(above[counts.indptr])


    def plot_sweep(sweep, chosen=None, title="", subtitle=""):
        """Cells keeping a guide (left) and guides per cell (right) against the threshold, a line per replicate."""
        colours = ["#006BA2", "#AA7C00", "#3EBCD2", "#A81829"]
        gridline = "#B7C6CF"
        fig, (left, right) = plt.subplots(1, 2, figsize=(7.5, 3.8), gridspec_kw={"wspace": 0.35})
        fig.subplots_adjust(left=0.09, right=0.88, bottom=0.14, top=0.68)
        panels = ((left, "with_guide", 100, "Cells with ≥ 1 guide, %"), (right, "moi", 1, "Guides per cell, mean"))
        for (name, rows), colour in zip(sweep.groupby("replicate", observed=True), colours):
            for ax, column, scale, _ in panels:
                ax.plot(rows["threshold"], scale * rows[column], color=colour, linewidth=1.5,
                        label=f"Replicate {name}")
            right.annotate(
                f"Rep {name}", (rows["threshold"].iloc[-1], rows["moi"].iloc[-1]),
                xytext=(4, 0), textcoords="offset points", va="center", fontsize=8, color=colour,
            )
        fig.legend(*left.get_legend_handles_labels(), loc="upper left", bbox_to_anchor=(0.08, 0.86),
                   ncol=4, frameon=False, fontsize=8.5, handlelength=1.5)
        for ax, _, _, label in panels:
            ax.set_title(label, loc="left", fontsize=9, color="0.3")
            ax.set_xlabel("Guide UMI threshold")
            ax.grid(axis="y", linewidth=0.5, color=gridline)
            ax.grid(axis="x", visible=False)
            ax.spines[["top", "right", "left"]].set_visible(False)
            ax.tick_params(axis="y", length=0)
            if chosen is not None:
                ax.axvline(chosen, color="0.1", linewidth=1, linestyle=(0, (3, 1)))
        fig.text(0.09, 0.94, title, fontsize=13, fontweight="bold", ha="left")
        fig.text(0.09, 0.885, subtitle, fontsize=9.5, color="0.3", ha="left")
        return fig

    return guides_at, plot_sweep


@app.cell
def cells_per_plot_fn(matplotlib, np, plt):
    def plot_cells_per(cells, unit="guide", bin_width=2, cap=None, title="", subtitle=""):
        """Histogram of cells per guide or per element, with the median; values past `cap` share the last bar."""
        blue, grey, gridline = "#006BA2", "#758D99", "#B7C6CF"
        cells = np.asarray(cells)
        top = cap if cap is not None else cells.max()
        edges = bin_width * np.arange(top // bin_width + 2)
        fig, ax = plt.subplots(figsize=(7.5, 3.6))
        fig.subplots_adjust(left=0.11, right=0.95, bottom=0.15, top=0.76)
        heights, _, bars = ax.hist(
            np.minimum(cells, top), bins=edges, color=blue, linewidth=0
        )
        ax.set(xlim=(0, edges[-1]), ylim=(0, heights.max() * 1.2))
        ax.set(xlabel=f"Cells assigned to the {unit}", ylabel=f"{unit.capitalize()}s")
        if cap is not None and (over := int((cells >= cap).sum())):
            bars[-1].set_color(grey)
            ax.annotate(
                f"{over:,} at {cap:,}+", (edges[-1], heights[-1]), xytext=(0, 3),
                textcoords="offset points", ha="right", fontsize=8, color="0.3",
            )
        median = np.median(cells)
        ax.axvline(median, ymax=1 / 1.2 + 0.08, color="0.1", linewidth=1)
        ax.annotate(
            f"median {median:,.0f} cells", (median, 1), xycoords=("data", "axes fraction"),
            xytext=(4, -9), textcoords="offset points", fontsize=8, color="0.1",
        )
        for axis in (ax.xaxis, ax.yaxis):
            axis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
        ax.grid(axis="y", linewidth=0.5, color=gridline)
        ax.grid(axis="x", visible=False)
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.tick_params(axis="y", length=0)
        fig.text(0.11, 0.93, title, fontsize=13, fontweight="bold", ha="left")
        fig.text(0.11, 0.875, subtitle, fontsize=9.5, color="0.3", ha="left")
        return fig


    return (plot_cells_per,)


@app.cell
def pairs_fn(np):
    def pair_table(pairs, pair_type, elements, readout_genes):
        """Label element-gene pairs and add the gene's symbol and the element's distance to its TSS.

        `pairs` holds `grna_target` (an `elements` id) and `response_id` (a readout gene id).
        The distance runs from the element's midpoint to the TSS, in the gene's orientation:
        negative is upstream.
        """
        out = (
            pairs[["grna_target", "response_id"]]
            .join(elements[["start", "end", "genomic_element"]], on="grna_target")
            .join(readout_genes[["symbol", "gene_tss", "gene_strand"]], on="response_id")
        )
        out["distance_to_tss"] = np.where(out["gene_strand"] == "-", -1, 1) * (
            (out["start"] + out["end"]) // 2 - out["gene_tss"]
        )
        return (
            out.assign(pair_type=pair_type)
            .rename(columns={"symbol": "gene_symbol"})
            [["grna_target", "response_id", "gene_symbol", "pair_type", "genomic_element", "distance_to_tss"]]
            .reset_index(drop=True)
        )


    return (pair_table,)


@app.cell
def guide_qc_fns(matplotlib, np, pd, plt):
    def house_style(ax, gridline="#B7C6CF"):
        """Value gridlines only, no top/right/left spines."""
        ax.grid(axis="y", linewidth=0.5, color=gridline)
        ax.grid(axis="x", visible=False)
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.tick_params(axis="y", which="both", length=0)


    def headline(fig, title, subtitle, left=0.11, top=0.93):
        """Title (the finding) and subtitle (what is measured), left-aligned above the axes."""
        fig.text(left, top, title, fontsize=13, fontweight="bold", ha="left")
        fig.text(left, top - 0.055, subtitle, fontsize=9.5, color="0.3", ha="left")


    def plot_depth(guides, gene_umis, cap=16, title="", subtitle=""):
        """Median gene UMIs by guides assigned, relative to one-guide cells, against a doublet's 2."""
        blue = "#006BA2"
        frame = pd.DataFrame({"guides": np.minimum(guides, cap), "gene": gene_umis})
        median = frame.groupby("guides")["gene"].median()
        ratio = median / median.loc[1]
        fig, ax = plt.subplots(figsize=(7.5, 3.6))
        fig.subplots_adjust(left=0.11, right=0.93, bottom=0.15, top=0.76)
        ax.axhline(2, color="0.1", linewidth=1, linestyle=(0, (4, 2)))
        ax.annotate("a doublet: twice one cell's UMIs", (ratio.index[0], 2), xytext=(0, 4),
                    textcoords="offset points", fontsize=8, color="0.1")
        ax.plot(ratio.index, ratio.to_numpy(), color=blue, linewidth=2, marker="o", markersize=4)
        ax.annotate(f"{ratio.iloc[-1]:.2f}×", (ratio.index[-1], ratio.iloc[-1]), xytext=(6, 0),
                    textcoords="offset points", va="center", fontsize=9, color=blue, fontweight="bold")
        ax.set_xticks(ratio.index, [str(v) for v in ratio.index[:-1]] + [f"{cap}+"])
        ax.set(xlim=(ratio.index[0] - 0.5, cap + 0.5), ylim=(0, 2.3), xlabel="Guides assigned to the cell",
               ylabel="Gene UMIs, × one-guide cells")
        ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
        ax.tick_params(axis="x", length=0, labelsize=8)
        house_style(ax)
        headline(fig, title, subtitle)
        return fig


    def calls_by_threshold(umis, top=200):
        """Guide calls (counts at or above the threshold) for every threshold from 1 to `top`."""
        tally = np.bincount(np.minimum(np.asarray(umis, dtype=np.int64), top + 1), minlength=top + 2)
        at_least = tally[::-1].cumsum()[::-1]
        return pd.Series(at_least[1 : top + 1], index=pd.RangeIndex(1, top + 1, name="threshold"), name="calls")


    def plot_calls_by_threshold(calls, chosen, title="", subtitle=""):
        """Guide calls against the threshold, from 2 UMIs, with the chosen threshold and the stretch from half to twice it."""
        blue = "#006BA2"
        calls = calls.loc[2:]
        fig, ax = plt.subplots(figsize=(7.5, 3.6))
        fig.subplots_adjust(left=0.11, right=0.95, bottom=0.15, top=0.76)
        ax.axvspan(chosen / 2, chosen * 2, color="#E9EDF0", linewidth=0)
        ax.plot(calls.index, calls.to_numpy() / 1e6, color=blue, linewidth=2)
        ax.axvline(chosen, color="0.1", linewidth=1)
        ax.annotate(f"{chosen:,.0f} UMIs: {calls.loc[int(chosen)] / 1e6:.2f}M calls", (chosen, 1),
                    xycoords=("data", "axes fraction"), xytext=(4, -9), textcoords="offset points",
                    fontsize=8, color="0.1")
        ticks = [t for t in (2, 5, 10, 20, 50, 100, 200, 500) if t <= calls.index[-1]]
        ax.set_xscale("log")
        ax.set_xticks(ticks, [str(t) for t in ticks])
        ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        ax.set(xlim=(2, calls.index[-1]), ylim=(0, calls.max() / 1e6 * 1.1),
               xlabel="UMI threshold", ylabel="Guide calls, millions")
        house_style(ax)
        headline(fig, title, subtitle)
        return fig


    def plot_by_group(values, groups, bin_width=2, unit="", xlabel="", title="", subtitle=""):
        """One histogram row per group, each as % of its own members, with its median."""
        blue = "#006BA2"
        values, groups = np.asarray(values), np.asarray(groups)
        names = sorted(set(groups))
        gap, row_height = 0.3, 0.75
        height = 1.5 + (row_height + gap) * len(names)
        bottom, top = 0.55 / height, 1 - 0.95 / height
        axes_height = (height * (top - bottom) - gap * (len(names) - 1)) / len(names)
        fig, axes = plt.subplots(len(names), 1, figsize=(7.5, height), sharex=True, sharey=True)
        fig.subplots_adjust(left=0.11, right=0.95, bottom=bottom, top=top, hspace=gap / axes_height)
        edges = bin_width * np.arange(values.max() // bin_width + 2)
        tallest = 0
        for ax, name in zip(axes, names):
            part = values[groups == name]
            heights, _, _ = ax.hist(part, bins=edges, weights=np.full(len(part), 100 / len(part)), color=blue, linewidth=0)
            tallest = max(tallest, heights.max())
            median = np.median(part)
            ax.axvline(median, color="0.1", linewidth=1)
            ax.annotate(f"median {median:,.0f}", (median, 1), xycoords=("data", "axes fraction"),
                        xytext=(4, -9), textcoords="offset points", fontsize=8, color="0.1")
            ax.set_title(f"{name} · {len(part):,} {unit}", loc="left", fontsize=9, color="0.3")
            house_style(ax)
        axes[0].set(xlim=(0, edges[-1]), ylim=(0, tallest * 1.3))
        axes[-1].set_xlabel(xlabel)
        fig.text(0.02, 0.5, "% of the row", rotation=90, va="center", fontsize=10)
        fig.text(0.11, 1 - 0.3 / height, title, fontsize=13, fontweight="bold", ha="left")
        fig.text(0.11, 1 - 0.52 / height, subtitle, fontsize=9.5, color="0.3", ha="left")
        return fig


    def plot_replicates(per_replicate, pairs, title="", subtitle=""):
        """Each guide's UMIs in one replicate against another, log axes, one panel per pair; one dot per guide."""
        blue = "#006BA2"
        short = matplotlib.ticker.FuncFormatter(
            lambda v, _: f"{v / 1e6:,.0f}M" if v >= 1e6 else f"{v / 1e3:,.0f}k" if v >= 1e3 else f"{v:,.0f}"
        )
        top = 10 ** np.ceil(np.log10(per_replicate.to_numpy().max()))
        fig, axes = plt.subplots(1, len(pairs), figsize=(7.5, 4.4))
        fig.subplots_adjust(left=0.11, right=0.95, bottom=0.13, top=0.78, wspace=0.35)
        for ax, (a, b) in zip(np.atleast_1d(axes), pairs):
            x, y = per_replicate[a].to_numpy(), per_replicate[b].to_numpy()
            both = (x > 0) & (y > 0)
            ax.scatter(x[both], y[both], s=2, color=blue, alpha=0.1, linewidths=0, rasterized=True)
            r = np.corrcoef(np.log10(x[both]), np.log10(y[both]))[0, 1]
            ax.set_title(f"r = {r:.2f} (log UMIs)", loc="left", fontsize=9, color="0.3")
            ax.set(xscale="log", yscale="log", xlim=(1, top), ylim=(1, top), aspect="equal",
                   xlabel=f"UMIs per guide, replicate {a}", ylabel=f"UMIs per guide, replicate {b}")
            for axis in (ax.xaxis, ax.yaxis):
                axis.set_major_locator(matplotlib.ticker.LogLocator(base=100))
                axis.set_major_formatter(short)
                axis.set_minor_locator(matplotlib.ticker.NullLocator())
            ax.spines[["top", "right"]].set_visible(False)
        fig.text(0.11, 0.94, title, fontsize=13, fontweight="bold", ha="left")
        fig.text(0.11, 0.895, subtitle, fontsize=9.5, color="0.3", ha="left")
        return fig


    def plot_one_channel(observed, expected, title="", subtitle=""):
        """Share of identical-guide groups with every cell in one channel, observed and expected from channel sizes."""
        blue, grey = "#006BA2", "#B7C6CF"
        fig, ax = plt.subplots(figsize=(7.5, 2.6))
        fig.subplots_adjust(left=0.31, right=0.93, bottom=0.2, top=0.66)
        values = [100 * expected, 100 * observed]
        ax.barh(["Expected from channel sizes", "Observed"], values, color=[grey, blue], height=0.6)
        for y, v in enumerate(values):
            ax.annotate(f"{v:.0f}%", (v, y), xytext=(4, 0), textcoords="offset points", va="center", fontsize=10)
        ax.set(xlim=(0, 100), xlabel="Groups with all their cells in one channel, %")
        ax.grid(axis="x", linewidth=0.5, color="#B7C6CF")
        ax.grid(axis="y", visible=False)
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.tick_params(axis="y", length=0)
        headline(fig, title, subtitle, left=0.02, top=0.9)
        return fig


    return (
        calls_by_threshold,
        headline,
        house_style,
        plot_by_group,
        plot_calls_by_threshold,
        plot_depth,
        plot_one_channel,
        plot_replicates,
    )


@app.cell
def pairwise_qc_plot_fn(headline, house_style, matplotlib, plt):
    def plot_pairwise_qc(pairs, trt_cut, cntrl_cut, title="", subtitle=""):
        """Each pair's treated against control cells with a nonzero count of its gene, log scales, the cuts dashed."""
        fig, ax = plt.subplots(figsize=(7.5, 5.4))
        fig.subplots_adjust(left=0.13, right=0.95, bottom=0.11, top=0.78)
        for passed, colour, label, alpha in ((True, "#006BA2", "pass", 0.25), (False, "#A81829", "fail", 0.9)):
            sub = pairs[pairs["pass_qc"] == passed]
            ax.scatter(
                sub["n_nonzero_trt"], sub["n_nonzero_cntrl"], s=5, color=colour, alpha=alpha,
                linewidths=0, rasterized=True, clip_on=False, zorder=3, label=f"{label}: {len(sub):,} pairs",
            )
        for line in (ax.axvline(trt_cut), ax.axhline(cntrl_cut)):
            line.set(color="0.1", linewidth=1, linestyle=(0, (4, 2)))
        for set_scale in (ax.set_xscale, ax.set_yscale):
            set_scale("symlog", linthresh=1, linscale=0.3)
        ax.set_xlim(0, pairs["n_nonzero_trt"].max() * 1.5)
        ax.set_ylim(0, pairs["n_nonzero_cntrl"].max() * 1.5)
        count = matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:,.0f}")
        ax.xaxis.set_major_formatter(count)
        ax.yaxis.set_major_formatter(count)
        ax.set(xlabel="Treated cells with a nonzero count", ylabel="Control cells with a nonzero count")
        ax.tick_params(labelsize=8)
        house_style(ax)
        ax.spines["left"].set_visible(True)
        legend = ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2, frameon=False, fontsize=9,
                            markerscale=2.5, handletextpad=0.2, borderaxespad=0.3)
        for handle in legend.legend_handles:
            handle.set_alpha(1)
        headline(fig, title, subtitle, left=0.13)
        return fig

    return (plot_pairwise_qc,)


if __name__ == "__main__":
    app.run()
