# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "marimo>=0.25.2",
#     "matplotlib==3.11.2",
#     "mudata==0.4.1",
#     "numpy==2.5.5",
#     "pandas==3.1.1",
#     "pysceptre[fast,io]==0.6.0",
#     "scipy==1.18.1",
# ]
# ///

import marimo

__generated_with = "0.25.2"
app = marimo.App(width="medium", auto_download=["html"])


@app.cell(hide_code=True)
def title(mo):
    mo.md(r"""
    # Analysing a CRISPR screen with pysceptre

    We test which candidate elements regulate which genes on day 0 of the
    hiPSC-to-endothelial-cell differentiation screen, read out on a targeted panel
    of 292 genes.
    """)
    return


@app.cell
def import_marimo():
    import marimo as mo

    return (mo,)


@app.cell
def imports(mo):
    import hashlib
    import urllib.request

    import matplotlib
    import matplotlib.pyplot as plt
    import mudata
    import numpy as np
    import pandas as pd
    import pysceptre
    from scipy import sparse, stats

    # marimo shows figures as SVG, with labels as real text: sharp at any zoom
    plt.rcParams["savefig.format"] = "svg"
    plt.rcParams["svg.fonttype"] = "none"
    plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial", "Helvetica"]

    mo.md(f"pysceptre **{pysceptre.__version__}**, mudata **{mudata.__version__}**")
    return (
        hashlib,
        matplotlib,
        mudata,
        np,
        pd,
        plt,
        pysceptre,
        sparse,
        stats,
        urllib,
    )


@app.cell
def settings(mo):
    # Cells to analyse, a random sample so a training run finishes in minutes.
    # None analyses every cell; only then do the numbers match the full screen.
    N_CELLS = 50_000
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

    The day-0 cells that passed QC, with their gene and guide counts, in one
    MuData file (`.h5mu`). The next cell downloads it once into `data/` and checks
    its md5.
    """)
    return


@app.cell(hide_code=True)
def get_data(hashlib, mo, urllib):
    DATA_URL = "https://mitra.stanford.edu/engreitz/oak/public/emattei/WTC11/day0_pass_qc.h5mu"
    DATA_MD5 = "0916e2f6f540fa6f3a83fa240a850102"
    data_path = mo.notebook_dir() / "data" / "day0_pass_qc.h5mu"

    if not data_path.exists():
        data_path.parent.mkdir(parents=True, exist_ok=True)
        _part = data_path.with_suffix(".part")
        _md5 = hashlib.md5()
        with urllib.request.urlopen(DATA_URL) as _response, _part.open("wb") as _out:
            with mo.status.progress_bar(
                total=int(_response.headers["Content-Length"]),
                title=f"Downloading {data_path.name}",
                show_rate=True,
                show_eta=True,
            ) as _bar:
                while _chunk := _response.read(8 << 20):
                    _out.write(_chunk)
                    _md5.update(_chunk)
                    _bar.update(increment=len(_chunk))
        if _md5.hexdigest() != DATA_MD5:
            _part.unlink()
            raise ValueError(f"md5 {_md5.hexdigest()}, expected {DATA_MD5}: the download is corrupt")
        _part.rename(data_path)

    mo.md(f"`{data_path}`: {data_path.stat().st_size / 1e9:.2f} GB")
    return (data_path,)


@app.cell(hide_code=True)
def md_load(mo):
    mo.md(r"""
    ## Load the screen

    Read the file with `mudata`. It holds two modalities over the same cells:
    `gene` (UMIs per gene) and `guide` (UMIs per guide, plus a `guide_assignment`
    layer saying which guides each cell carries).
    """)
    return


@app.cell
def load(N_CELLS, SAMPLE_SEED, data_path, mo, mudata, np):
    mdata = mudata.read_h5mu(data_path)
    _n_all = mdata.n_obs
    if N_CELLS is not None and N_CELLS < _n_all:
        # sorted, so the sampled cells keep their order in the file; the same cells every run
        _keep = np.sort(np.random.default_rng(SAMPLE_SEED).choice(_n_all, N_CELLS, replace=False))
        mdata = mdata[_keep].copy()

    mo.vstack(
        [
            mo.md(
                f"**{mdata.n_obs:,}** of the {_n_all:,} cells"
                + (", a random sample (`N_CELLS`)" if mdata.n_obs < _n_all else "")
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

    The file keeps only the cells that passed QC, so the whole curves come from
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
        title="The 1,000-UMI cut falls on every channel's cliff",
        subtitle="Day 0: barcodes ranked by gene UMIs, log scales; dashed line, the cut",
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
                    f"Day 0, {len(readout_share):,} cells: share of each cell's gene UMIs "
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
            f"Day 0: share of each cell's gene UMIs from the {int(genes['readout'].sum())} "
            f"readout genes; {int((readout_share < 0.5).sum())} cells below 50% not shown"
        ),
    )
    return


@app.cell(hide_code=True)
def md_readout_by_channel(mo):
    mo.md(r"""
    ### By channel

    Pick a batch to see its channels, one row each. In batch 4, the second hump
    is a single channel.
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
        subtitle="Day 0: share of each cell's gene UMIs from the readout genes, one row per channel",
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
        subtitle="Day 0: barcodes ranked by guide UMIs, log scales",
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
        subtitle="Day 0: cells keeping a guide, and guides per cell, by guide UMI threshold",
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
            f"Day 0, {len(guides_per_cell):,} cells: guides assigned per cell "
            f"({assignment.threshold:.0f}-UMI threshold); mean {guides_per_cell.mean():.1f}"
        ),
    )
    return (guides_per_cell,)


@app.cell(hide_code=True)
def note_low_moi(mo):
    mo.md(r"""
    /// note | Later: an internal low-MOI control

    About a quarter of the cells carry exactly one called guide. Analysed on their
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
            "Day 0: median gene UMIs per cell, relative to cells with one guide "
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
                subtitle="Day 0: guide calls, a guide at or above the threshold in a cell; shaded, half to twice it",
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
    chosen above, all 36 channels pooled.
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
            f"Day 0: cells assigned to each of the {len(cells_per_guide):,} guides "
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
            f"Day 0: cells per guide by guide type, % of the type's guides "
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
                subtitle="Day 0: total UMIs of each guide in each replicate, one dot per guide",
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
                subtitle=f"Day 0: {len(identical_sets):,} groups, each of cells carrying exactly the same 3 or more guides",
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
            f"Day 0: cells carrying any guide of each of the {len(elements):,} targeted elements "
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
    return MAX_DISTANCE, discovery_pairs


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
    sides: at least 7 **treated** cells (carrying any of the element's guides) and
    7 **control** cells (all the other cells) with a nonzero count, sceptre's
    defaults. pysceptre's power check applies this itself; discovery expects pairs
    that already pass, so we apply it here, with the same function, to every pair.
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
    N_NONZERO_TRT, N_NONZERO_CNTRL = 7, 7  # sceptre's defaults

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
                "all pass: day 0 barely expresses them."
            ),
            mo.ui.table(_by_gene.loc[_failing].sort_values("fail_qc", ascending=False).round(2)),
        ]
    )
    return N_NONZERO_CNTRL, N_NONZERO_TRT, pairs_qc


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
    return (covariate_matrix,)


@app.cell
def md_calibration(mo):
    mo.md(r"""
    ## Calibration check

    Before testing real pairs, we ask whether the test invents effects. pysceptre
    builds fake targets out of control guides, as many guides per fake target as a
    median element has, pairs them with readout genes, and runs exactly the test the
    discovery run will use. No fake target is real, so the p-values should be
    uniform between 0 and 1, and almost nothing should pass the multiple-testing
    correction. Too many small p-values mean the test, with these covariates and
    settings, would call false hits in discovery too.

    The check doesn't set anything for discovery: it is only meaningful if both runs
    get the same settings, so they are defined once here. We run it twice, with two
    kinds of control guide:

    - **non-targeting**: guides with no site in the genome;
    - **safe-targeting**: guides that cut or bind the genome, but away from every
      element. A safe-targeting guide close enough to a gene to repress it would show
      up as a real effect, not as noise.
    """)
    return


@app.cell
def calibration_settings(
    assignment,
    elements,
    guide_cells,
    guides,
    mo,
    np,
    pairs_qc,
):
    SIDE = "both"  # two-sided: decreases and increases
    SEED = 1
    FDR = 0.2  # Benjamini-Hochberg level, the DC-TAP-seq paper's, for the calibration and for discovery
    N_JOBS = 8  # the box's CPUs

    # R's defaults: as many pairs as discovery will test, and the median element's number of guides
    _tested = pairs_qc["pair_type"] == "discovery"
    N_CALIBRATION_PAIRS = int((pairs_qc["pass_qc"] & _tested).sum())
    PASS_QC_RATE = N_CALIBRATION_PAIRS / int(_tested.sum())
    CALIBRATION_GROUP_SIZE = int(elements["guides"].median())

    control_guide_cells = {
        _kind: guide_cells(_kind, guides, assignment.assigned) for _kind in ("non-targeting", "safe-targeting")
    }

    mo.md(
        f"- `SIDE` = `{SIDE}`, `moi=\"high\"`: each fake target against every other cell\n"
        f"- `FDR` = {FDR}: the Benjamini-Hochberg level, here and in discovery\n"
        f"- **{N_CALIBRATION_PAIRS:,}** calibration pairs, as many as the discovery pairs that pass QC "
        f"({PASS_QC_RATE:.1%} of them)\n"
        f"- **{CALIBRATION_GROUP_SIZE}** guides per fake target, the median element's\n"
        + "\n".join(
            f"- {_kind}: {len(_cells):,} guides, a median of "
            f"{np.median([len(_c) for _c in _cells.values()]):.0f} cells each"
            for _kind, _cells in control_guide_cells.items()
        )
    )
    return (
        CALIBRATION_GROUP_SIZE,
        FDR,
        N_CALIBRATION_PAIRS,
        N_JOBS,
        PASS_QC_RATE,
        SEED,
        SIDE,
        control_guide_cells,
    )


@app.cell
def calibration_nt(
    CALIBRATION_GROUP_SIZE,
    FDR,
    N_CALIBRATION_PAIRS,
    N_JOBS,
    N_NONZERO_CNTRL,
    N_NONZERO_TRT,
    PASS_QC_RATE,
    SEED,
    SIDE,
    control_guide_cells,
    covariate_matrix,
    gene_ids,
    mo,
    pysceptre,
    response_matrix,
    stats,
):
    calibration_nt = pysceptre.run_calibration_check(
        response_matrix,
        gene_ids,
        covariate_matrix,
        control_guide_cells["non-targeting"],
        n_calibration_pairs=N_CALIBRATION_PAIRS,
        calibration_group_size=CALIBRATION_GROUP_SIZE,
        n_nonzero_trt_thresh=N_NONZERO_TRT,
        n_nonzero_cntrl_thresh=N_NONZERO_CNTRL,
        pass_qc_rate=PASS_QC_RATE,
        side=SIDE,
        moi="high",
        seed=SEED,
        n_jobs=N_JOBS,
    )
    _q = stats.false_discovery_control(calibration_nt["p_value"], method="bh")
    mo.md(
        f"**non-targeting**: {len(calibration_nt):,} pairs tested, "
        f"{calibration_nt['grna_target'].nunique():,} fake targets. "
        f"{(calibration_nt['p_value'] < 0.05).mean():.1%} have p < 0.05 (5% expected); "
        f"**{int((_q < FDR).sum())}** pass Benjamini-Hochberg at {FDR} (0 expected)."
    )
    return (calibration_nt,)


@app.cell
def calibration_safe(
    CALIBRATION_GROUP_SIZE,
    FDR,
    N_CALIBRATION_PAIRS,
    N_JOBS,
    N_NONZERO_CNTRL,
    N_NONZERO_TRT,
    PASS_QC_RATE,
    SEED,
    SIDE,
    control_guide_cells,
    covariate_matrix,
    gene_ids,
    mo,
    pysceptre,
    response_matrix,
    stats,
):
    calibration_safe = pysceptre.run_calibration_check(
        response_matrix,
        gene_ids,
        covariate_matrix,
        control_guide_cells["safe-targeting"],
        n_calibration_pairs=N_CALIBRATION_PAIRS,
        calibration_group_size=CALIBRATION_GROUP_SIZE,
        n_nonzero_trt_thresh=N_NONZERO_TRT,
        n_nonzero_cntrl_thresh=N_NONZERO_CNTRL,
        pass_qc_rate=PASS_QC_RATE,
        side=SIDE,
        moi="high",
        seed=SEED,
        n_jobs=N_JOBS,
    )
    _q = stats.false_discovery_control(calibration_safe["p_value"], method="bh")
    mo.md(
        f"**safe-targeting**: {len(calibration_safe):,} pairs tested, "
        f"{calibration_safe['grna_target'].nunique():,} fake targets. "
        f"{(calibration_safe['p_value'] < 0.05).mean():.1%} have p < 0.05 (5% expected); "
        f"**{int((_q < FDR).sum())}** pass Benjamini-Hochberg at {FDR} (0 expected)."
    )
    return


@app.cell
def md_calibration_hits(mo):
    mo.md(r"""
    ### Which guide makes a non-targeting hit?

    A fake target is 15 non-targeting guides pooled. If one of them really acts on
    the gene it is paired with, the hit is a real effect, not a failure of the test.
    For each hit, the gene's expression in each guide's cells, relative to all
    cells (UMIs scaled by each cell's gene UMIs): a guide near 1 does nothing, a
    guide far from 1 carries the effect.
    """)
    return


@app.cell
def calibration_hits(
    FDR,
    calibration_nt,
    cell_ids,
    control_guide_cells,
    gene_ids,
    mdata,
    mo,
    np,
    pd,
    readout_genes,
    response_matrix,
    stats,
):
    _q = stats.false_discovery_control(calibration_nt["p_value"], method="bh")
    _hits = calibration_nt.loc[_q < FDR].sort_values("p_value")
    _depth = mdata["gene"].obs.loc[cell_ids, "total_counts"].to_numpy()
    _depth = _depth / np.median(_depth)

    _rows = []
    for _pair in _hits.itertuples():
        _expression = response_matrix[gene_ids.index(_pair.response_id)].toarray().ravel() / _depth
        for _guide in _pair.grna_target.split("&"):
            _cells = control_guide_cells["non-targeting"][_guide]
            _rows.append(
                {
                    "gene": readout_genes.loc[_pair.response_id, "symbol"],
                    "pair_fold_change": _pair.fold_change,
                    "p_value": _pair.p_value,
                    "guide": _guide,
                    "guide_cells": len(_cells),
                    "guide_fold_change": _expression[_cells].mean() / _expression.mean(),
                }
            )
    calibration_hit_guides = pd.DataFrame(
        _rows, columns=["gene", "pair_fold_change", "p_value", "guide", "guide_cells", "guide_fold_change"]
    )
    # per hit, the guide furthest from 1 on a log scale
    _far = np.abs(np.log2(calibration_hit_guides["guide_fold_change"].clip(lower=1e-3)))
    _drivers = calibration_hit_guides.loc[_far.groupby([calibration_hit_guides["gene"], calibration_hit_guides["p_value"]]).idxmax()]

    mo.vstack(
        [
            mo.md(f"**{len(_hits)}** non-targeting pairs pass Benjamini-Hochberg at {FDR}. The guide that moves the gene most, per hit:"),
            _drivers.sort_values("p_value").reset_index(drop=True),
            mo.ui.table(calibration_hit_guides, selection=None),
        ]
        if len(_hits)
        else [mo.md(f"**0** non-targeting pairs pass Benjamini-Hochberg at {FDR}: no hit to trace to a guide.")]
    )
    return


@app.cell
def md_power_check(mo):
    mo.md(r"""
    ## Power check

    The calibration asked whether the test invents effects; the power check asks
    whether it finds the ones it should. It tests two kinds of positive control,
    each expected to have a negative effect on its gene, and reports them apart:

    - the pairs the design names: promoters and validated enhancers against their
      genes;
    - the host-gene pairs: each element inside a gene's body against that gene.

    The settings are the calibration's, and will be discovery's. Pairs that fail QC
    are reported, not dropped: dropping them would hide the controls the screen had
    too few cells to test.
    """)
    return


@app.cell
def power_check(
    N_JOBS,
    N_NONZERO_CNTRL,
    N_NONZERO_TRT,
    SEED,
    SIDE,
    covariate_matrix,
    gene_ids,
    grna_target_cells,
    host_gene_pairs,
    mo,
    pd,
    positive_controls,
    pysceptre,
    response_matrix,
):
    # the design's positive controls, and the host-gene pairs
    _controls = pd.concat([positive_controls, host_gene_pairs], ignore_index=True)

    power_check = pysceptre.run_power_check(
        response_matrix,
        gene_ids,
        covariate_matrix,
        grna_target_cells,
        positive_control_pairs=_controls[["response_id", "grna_target"]],
        n_nonzero_trt_thresh=N_NONZERO_TRT,
        n_nonzero_cntrl_thresh=N_NONZERO_CNTRL,
        side=SIDE,
        moi="high",
        seed=SEED,
        n_jobs=N_JOBS,
    ).merge(
        _controls[["grna_target", "response_id", "gene_symbol", "pair_type", "genomic_element"]],
        on=["grna_target", "response_id"],
    )
    _tested = power_check[power_check["pass_qc"]]
    _knockdown = (_tested["p_value"] < 0.05) & (_tested["fold_change"] < 1)

    mo.vstack(
        [
            mo.md(
                f"**{len(_tested):,}** of {len(power_check):,} pairs pass QC and are tested. "
                + "; ".join(
                    f"{_type.replace('_', '-')} pairs: **{int(_k.sum()):,}** of {len(_k):,} ({_k.mean():.0%}) "
                    "have a negative effect at p < 0.05"
                    for _type, _k in _knockdown.groupby(_tested["pair_type"], sort=False)
                )
                + ". No multiple-testing correction, as in sceptre: these measure power, they are not discoveries."
            ),
            _tested.assign(negative_effect=_knockdown.to_numpy())
            .groupby(["pair_type", "genomic_element"], observed=True, sort=False)
            .agg(tested=("p_value", "size"), negative_effect=("negative_effect", "sum"), median_fold_change=("fold_change", "median")),
            mo.ui.table(power_check, selection=None),
        ]
    )
    return


@app.cell
def md_discovery_run(mo):
    mo.md(r"""
    ## Discovery

    The discovery pairs that passed QC, tested with the calibration's settings.
    Positive controls and host-gene pairs are the power check's.

    pysceptre returns a p-value per pair and corrects none of them, so the
    correction is ours, in the next cell.
    """)
    return


@app.cell
def discovery_run(
    N_JOBS,
    SEED,
    SIDE,
    covariate_matrix,
    discovery_pairs,
    gene_ids,
    grna_target_cells,
    pairs_qc,
    pysceptre,
    response_matrix,
):
    # discovery_pairs, keeping those that pass QC
    _passing = pairs_qc.loc[pairs_qc["pass_qc"], ["grna_target", "response_id"]]
    _tested = discovery_pairs[["response_id", "grna_target"]].merge(_passing, on=["grna_target", "response_id"])

    # the calibration's settings: the check vouches only for a run made with them
    discovery_results = pysceptre.run_discovery_analysis(
        response_matrix,
        gene_ids,
        covariate_matrix,
        grna_target_cells,
        _tested,
        side=SIDE,
        moi="high",
        seed=SEED,
        n_jobs=N_JOBS,
    )
    return (discovery_results,)


@app.cell
def md_discovery_table(mo):
    mo.md(r"""
    ### Hits

    Each tested pair with its annotations, and a q-value from Benjamini-Hochberg at
    `FDR` (0.2, the DC-TAP-seq paper's level) over the pairs with `in_fdr`. Self
    promoters keep their p-value but get no q-value. A hit is a pair with q below
    `FDR`; its fold change says whether the effect is negative (below 1) or positive
    (above 1).
    """)
    return


@app.cell
def discovery_table(FDR, discovery_results, mo, np, pairs_qc, stats):
    # pysceptre's results with each pair's annotations and QC counts
    discovery_table = discovery_results.merge(
        pairs_qc.drop(columns="pass_qc"), on=["grna_target", "response_id"], how="left", validate="one_to_one"
    )
    assert discovery_table["pair_type"].eq("discovery").all()  # only discovery pairs, each once
    assert discovery_table["p_value"].notna().all()

    # q-values only for the pairs in the correction: all but the self promoters
    _in_fdr = discovery_table["in_fdr"].fillna(False).to_numpy(dtype=bool)
    discovery_table["q_value"] = np.nan
    discovery_table.loc[_in_fdr, "q_value"] = stats.false_discovery_control(
        discovery_table.loc[_in_fdr, "p_value"], method="bh"
    )

    _hits = discovery_table.loc[discovery_table["q_value"] < FDR].sort_values("p_value")
    _down = int((_hits["fold_change"] < 1).sum())
    _summary = discovery_table.groupby("genomic_element", observed=True).agg(
        tested=("p_value", "size"),
        p_below_005=("p_value", lambda _p: int((_p < 0.05).sum())),
        in_fdr=("in_fdr", lambda _f: int(_f.fillna(False).sum())),
        hits=("q_value", lambda _q: int((_q < FDR).sum())),
    )

    mo.vstack(
        [
            mo.md(
                f"**{len(discovery_table):,} pairs tested.** Of the {int(_in_fdr.sum()):,} in the correction, "
                f"**{len(_hits)}** pass Benjamini-Hochberg at {FDR}: **{_down}** negative effects, "
                f"{len(_hits) - _down} positive. Outside the correction, `p_below_005` can be read "
                "against chance: 5% of the pairs tested."
            ),
            _summary,
            mo.ui.table(
                _hits[
                    [
                        "grna_target",
                        "gene_symbol",
                        "genomic_element",
                        "distance_to_tss",
                        "n_nonzero_trt",
                        "fold_change",
                        "pct_change_es_ci_low",
                        "pct_change_es_ci_high",
                        "p_value",
                        "q_value",
                    ]
                ],
                selection=None,
            ),
        ]
    )
    return (discovery_table,)


@app.cell
def md_pair_classes(mo):
    mo.md(r"""
    ### Hits by pair class and distance

    A pair's class is where its element sits: on a protein-coding TSS
    (`promoter-gene`), in an exon or an intron, or outside genes
    (`intergenic-gene`, the guide table's `distal element`). Validated enhancers
    keep their own class. The distance runs from the element's midpoint to the
    gene's TSS. Only the pairs in the correction are counted.
    """)
    return


@app.cell
def pair_classes(FDR, MAX_DISTANCE, discovery_table, mo, pd):
    PAIR_CLASSES = {
        "promoter": "promoter-gene",
        "exon": "exon-gene",
        "intron": "intron-gene",
        "distal element": "intergenic-gene",
        "enhancer": "enhancer-gene",
    }
    _in_fdr = discovery_table.loc[discovery_table["in_fdr"].fillna(False).to_numpy(dtype=bool)]
    _pairs = pd.DataFrame(
        {
            "pair_class": pd.Categorical(
                _in_fdr["genomic_element"].astype(str).map(PAIR_CLASSES), categories=list(PAIR_CLASSES.values())
            ),
            "distance": pd.cut(
                _in_fdr["distance_to_tss"].abs(),
                [0, 10_000, 100_000, MAX_DISTANCE],
                labels=["≤ 10 kb", "10–100 kb", "100 kb–1 Mb"],
                include_lowest=True,
            ),
            "distance_kb": _in_fdr["distance_to_tss"].abs() / 1000,
            "hit": (_in_fdr["q_value"] < FDR).to_numpy(),
            "down": ((_in_fdr["q_value"] < FDR) & (_in_fdr["fold_change"] < 1)).to_numpy(),
        }
    )
    pair_class_hits = _pairs[_pairs["hit"]]
    _negative = pair_class_hits[pair_class_hits["down"]]
    _positive = pair_class_hits[~pair_class_hits["down"]]

    pair_class_counts = _pairs.groupby("pair_class", observed=False).agg(
        tested=("hit", "size"),
        hits=("hit", "sum"),
        negative_effect=("down", "sum"),
    )
    pair_class_counts["positive_effect"] = pair_class_counts["hits"] - pair_class_counts["negative_effect"]
    pair_class_counts["median_hit_distance_kb"] = pair_class_hits.groupby("pair_class", observed=False)["distance_kb"].median().round(1)

    mo.vstack(
        [
            mo.md(f"**{int(_pairs['hit'].sum())}** hits at FDR {FDR}, by class:"),
            pair_class_counts,
            mo.md("Negative effects by class and distance to the TSS:"),
            pd.crosstab(_negative["pair_class"], _negative["distance"], margins=True, margins_name="all", dropna=False),
            mo.md("Positive effects by class and distance to the TSS:"),
            pd.crosstab(_positive["pair_class"], _positive["distance"], margins=True, margins_name="all", dropna=False),
        ]
    )
    return PAIR_CLASSES, pair_class_counts, pair_class_hits


@app.cell
def pair_classes_plot(
    FDR,
    pair_class_counts,
    pair_class_hits,
    plot_hits_by_class,
):
    _top = pair_class_counts["hits"].nlargest(2)

    plot_hits_by_class(
        pair_class_hits,
        tested=pair_class_counts["tested"],
        title=f"{_top.sum() / pair_class_counts['hits'].sum():.0%} of hits are {_top.index[0]} or {_top.index[1]} pairs",
        subtitle=f"Day 0: hits at FDR {FDR} by pair class, direction of effect and distance to the TSS",
    )
    return


@app.cell
def md_filter_sweep(mo):
    mo.md(r"""
    ### Fewer, better-powered pairs in the correction

    A pair whose gene is barely expressed in the element's cells cannot reach
    significance, but it still counts in the correction. Dropping such pairs raises
    the cutoff for the rest, as long as the filter does not depend on the result.
    Two filters, at the same thresholds, over the pairs with `in_fdr`:

    - **observed**: `n_nonzero_trt`, the treated cells that express the gene. It
      looks at the treated cells, so near the threshold it keeps the pairs whose
      expression rose by chance and drops those where it fell;
    - **expected**: the element's treated cells × the share of all cells that
      express the gene, which is what `n_nonzero_trt` would be if the element did
      nothing. It never looks at which cells are treated, so it is independent of
      the test when there is no effect.

    Set the threshold on power grounds, before reading the hits: picking the one
    that gives the most hits undoes the FDR control.
    """)
    return


@app.cell
def filter_sweep(
    FDR,
    N_NONZERO_TRT,
    discovery_table,
    gene_ids,
    grna_target_cells,
    mo,
    np,
    pd,
    response_matrix,
    stats,
):
    # what n_nonzero_trt would be if the element did nothing: its treated cells × the gene's detection rate
    _detection = pd.Series(np.diff(response_matrix.indptr) / response_matrix.shape[1], index=gene_ids)
    _treated = pd.Series({_target: len(_cells) for _target, _cells in grna_target_cells.items()})
    _in_fdr = discovery_table.loc[discovery_table["in_fdr"].fillna(False).to_numpy(dtype=bool)]
    _in_fdr = _in_fdr.assign(
        expected_nonzero_trt=_treated.loc[_in_fdr["grna_target"]].to_numpy()
        * _detection.loc[_in_fdr["response_id"]].to_numpy()
    )


    def _hits(_pairs):
        return int((stats.false_discovery_control(_pairs["p_value"], method="bh") < FDR).sum())


    _rows = []
    for _min_cells in [N_NONZERO_TRT, 10, 15, 20, 30, 50, 100, 200]:
        _observed = _in_fdr[_in_fdr["n_nonzero_trt"] >= _min_cells]
        _expected = _in_fdr[_in_fdr["expected_nonzero_trt"] >= _min_cells]
        _rows.append(
            {
                "min_cells": _min_cells,
                "pairs_observed": len(_observed),
                "hits_observed": _hits(_observed),
                "pairs_expected": len(_expected),
                "hits_expected": _hits(_expected),
            }
        )
    filter_sweep = pd.DataFrame(_rows).set_index("min_cells")

    mo.vstack(
        [
            mo.md(
                f"Hits at FDR {FDR} as pairs with fewer expressing treated cells leave the correction. "
                f"At {N_NONZERO_TRT}, the pairwise QC already applied, the observed filter keeps all "
                f"{filter_sweep.loc[N_NONZERO_TRT, 'pairs_observed']:,} pairs."
            ),
            filter_sweep,
        ]
    )
    return


@app.cell
def note_ihw(mo):
    mo.md(r"""
    /// note | Later: weight pairs by power instead of cutting them

    The sweep drops pairs below a threshold, and picking the threshold from its hits
    would undo the FDR control. Independent hypothesis weighting (IHW; Ignatiadis
    et al., 2016) uses the same covariate, the expected expressing treated cells, to
    give better-powered pairs more of the error budget: no threshold to choose, and
    the FDR still controlled. It is an R package; this needs a Python version.
    ///
    """)
    return


@app.cell
def md_power(mo):
    mo.md(r"""
    ## Power: which pairs could the screen detect?

    A pair that is not a hit either has no effect or could not have been detected.
    Power separates the two: if this element had a negative effect of a given size
    on this gene, how likely is this screen to call it? pysceptre answers it per
    pair in closed form (its port of PerturbPlan's post-hoc power), from three
    things the screen already holds:

    - each gene's expression and dispersion, from the same negative-binomial fits
      the test uses;
    - the cells each of the element's guides was assigned to;
    - the nominal p-value that Benjamini-Hochberg at `FDR` actually applied here.

    Effects are negative, at 5, 10, 15, 25 and 50%; the DC-TAP-seq paper reports power
    for a 15% effect. A pair counts as powered at 80%. The across-guide spread of
    the fold change, `FOLD_CHANGE_SD`, is not measured here: 0.13 is the value
    PerturbPlan's example uses and pysceptre was validated with.
    """)
    return


@app.cell
def gene_fits(
    N_JOBS,
    covariate_matrix,
    gene_ids,
    mo,
    pysceptre,
    response_matrix,
):
    # each gene's negative-binomial fit, as discovery made it: the expression scale the power formula needs
    gene_fits = pysceptre.pipeline.discovery.fit_all_genes(
        response_matrix, gene_ids, covariate_matrix, batch_width=1, n_jobs=N_JOBS
    )
    mo.md(f"{len(gene_fits)} genes fitted")
    return (gene_fits,)


@app.cell
def pair_power(
    FDR,
    PAIR_CLASSES,
    assignment,
    covariate_matrix,
    discovery_table,
    gene_fits,
    gene_ids,
    guides,
    mo,
    np,
    pd,
    pysceptre,
    response_matrix,
):
    EFFECT_SIZES = [0.05, 0.10, 0.15, 0.25, 0.50]  # negative effects; 15% is the DC-TAP-seq paper's reference
    FOLD_CHANGE_SD = 0.13  # across-guide sd of the fold change: PerturbPlan's example value, not measured here
    MIN_POWER = 0.8

    # one row per guide and its element, and the cells each guide was assigned to
    _targeted = np.flatnonzero(guides["intended_target_chr"].notna().to_numpy())
    _g = guides.iloc[_targeted]
    _design = pd.DataFrame(
        {
            "grna_id": _g.index,
            "grna_target": _g["intended_target_chr"].astype(str)
            + ":"
            + _g["intended_target_start"].astype("int64").astype(str)
            + "-"
            + _g["intended_target_end"].astype("int64").astype(str),
        }
    )
    _calls = assignment.assigned.tocsr()
    _cells_per_grna = pysceptre.analytical_power.cells_per_grna_from_assignments(
        _design, {guides.index[_i]: _calls.indices[_calls.indptr[_i] : _calls.indptr[_i + 1]] for _i in _targeted}
    )
    _baseline = pysceptre.analytical_power.baseline_expression_stats_from_fits(
        covariate_matrix, gene_fits, gene_subset=gene_ids
    )

    # the nominal p-value Benjamini-Hochberg at FDR applied; the left tail at half of it, as pysceptre validated
    _tested = discovery_table.loc[discovery_table["in_fdr"].fillna(False).to_numpy(dtype=bool)]
    POWER_CUTOFF = pysceptre.analytical_power.bh_nominal_cutoff(_tested["p_value"].to_numpy(), FDR)

    pair_power = _tested[["grna_target", "response_id", "gene_symbol", "genomic_element", "distance_to_tss", "p_value", "q_value"]].copy()
    for _es in EFFECT_SIZES:
        pair_power[f"power_{round(100 * _es)}pct"] = pysceptre.analytical_power.compute_power(
            _tested[["grna_target", "response_id"]],
            _cells_per_grna,
            _baseline,
            fold_change_mean=1 - _es,
            fold_change_sd=FOLD_CHANGE_SD,
            cutoff=POWER_CUTOFF / 2,
            num_total_cells=response_matrix.shape[1],
            side="left",
        )["power"].to_numpy()

    _columns = [f"power_{round(100 * _es)}pct" for _es in EFFECT_SIZES]
    _class = pd.Categorical(pair_power["genomic_element"].astype(str).map(PAIR_CLASSES), categories=list(PAIR_CLASSES.values()))
    _powered = (pair_power[_columns] >= MIN_POWER).groupby(_class, observed=False).mean().mul(100).round(1)
    _powered.insert(0, "pairs", pd.Series(_class).value_counts().reindex(_powered.index).to_numpy())
    _hit = (pair_power["q_value"] < FDR).to_numpy()
    _well = (pair_power["power_15pct"] >= MIN_POWER).to_numpy()

    mo.vstack(
        [
            mo.md(
                f"Benjamini-Hochberg at {FDR} applied p < {POWER_CUTOFF:.2e}. "
                f"**{_well.mean():.0%}** of the {len(pair_power):,} pairs in the correction have "
                f"{MIN_POWER:.0%} power for a 15% negative effect; "
                f"{int((_hit & _well).sum())} of the {int(_hit.sum())} hits are among them."
            ),
            mo.md(f"Share of pairs with {MIN_POWER:.0%} power, by class and effect size (%):"),
            _powered.rename_axis("pair class"),
        ]
    )
    return EFFECT_SIZES, MIN_POWER, pair_power


@app.cell(hide_code=True)
def md_power_curves(mo):
    mo.md(r"""
    ### How many pairs could the screen detect?

    Power depends most on how much the gene is expressed, so the pairs are split
    into five groups by their gene's mean UMIs per cell, a fifth of the pairs each.
    Each line counts a group's pairs with 80% power at each effect size. At no
    effect, power is the test's false-positive rate, so no pair is powered.
    """)
    return


@app.cell
def power_curves(
    EFFECT_SIZES,
    MIN_POWER,
    gene_ids,
    mo,
    np,
    pair_power,
    pd,
    plot_powered_by_effect,
    response_matrix,
):
    # genes in five groups of mean expression, cut so each group holds a fifth of the pairs
    _mean = pd.Series(np.asarray(response_matrix.mean(axis=1)).ravel(), index=list(gene_ids))
    _expression = pair_power["response_id"].map(_mean)
    _cuts = np.quantile(_expression, [0.2, 0.4, 0.6, 0.8])
    _names = [f"< {_cuts[0]:.2g}"] + [f"{_a:.2g}–{_b:.2g}" for _a, _b in zip(_cuts[:-1], _cuts[1:])] + [f"≥ {_cuts[-1]:.2g}"]
    _group = pd.cut(_expression, [-np.inf, *_cuts, np.inf], right=False, labels=_names)

    # rows: effect sizes; columns: expression groups; values: pairs with MIN_POWER
    _powered = pd.DataFrame(
        {
            _es: (pair_power[f"power_{round(100 * _es)}pct"] >= MIN_POWER).groupby(_group, observed=False).sum()
            for _es in EFFECT_SIZES
        }
    ).T
    _totals = pd.Series(_group).value_counts()
    _table = _powered.T.rename(columns=lambda _es: f"{_es:.0%}")
    _table.insert(0, "pairs", _totals.reindex(_table.index).to_numpy())

    mo.vstack(
        [
            plot_powered_by_effect(
                _powered,
                totals=_totals,
                reference=0.15,
                reference_label="the DC-TAP-seq paper's 15%",
                legend_title="Gene expression, mean UMIs per cell",
                title="Pairs of highly expressed genes are powered at smaller effects",
                subtitle=f"Day 0: pairs in the correction with {MIN_POWER:.0%} power, by size of negative effect",
            ),
            mo.md(f"Pairs with {MIN_POWER:.0%} power, by gene mean UMIs per cell and negative effect:"),
            _table.rename_axis("mean UMIs per cell"),
        ]
    )
    return


@app.cell(hide_code=True)
def md_long_range(mo):
    mo.md(r"""
    ## Long-range links: do they clear the trans level?

    An element cannot act directly on a gene on another chromosome, so its trans
    effects are indirect, through some gene it does regulate, or false calls. Their
    rate is what an element-gene pair reaches without direct regulation, and a cis
    link is credible only where the cis rate clears that trans level. Near the TSS
    it clears it by far; at long range, where direct links are rare, it may not.

    Every element in the correction is tested against the panel's genes on other
    chromosomes, with discovery's own test, fits and seed. pysceptre's
    `run_specificity_check` then compares the two rates gene for gene in each
    distance bin. The **trans/cis ratio** is the share of a bin's links the trans
    level alone would produce: near 0 the links are real, near 1 they are no more
    frequent than indirect effects.

    There are about twenty times as many trans pairs as cis pairs, so they run in
    batches. pysceptre fits the elements in chunks, and a trans test reproduces
    discovery only if each element is fit with the same chunk-mates. `trans_fit_chunks`
    rebuilds discovery's chunks from its pairs, and `fit_chunk_ids` gives each call
    a few whole chunks: a call tests only the elements in its chunks and refuses any
    other. Every batch returns exactly what one call over all the pairs returns, so
    the batches can run one after another, as here, or as separate jobs.

    Only negative effects are counted: the trans test does not resample a pair
    whose fold change is at least 1, so the positive hits have no trans level here.
    """)
    return


@app.cell
def trans_results(
    N_JOBS,
    SEED,
    SIDE,
    covariate_matrix,
    discovery_pairs,
    elements,
    gene_fits,
    gene_ids,
    genes,
    grna_target_cells,
    mo,
    np,
    pairs_qc,
    pd,
    pysceptre,
    response_matrix,
):
    TRANS_CHUNKS_PER_CALL = 4  # discovery's fit chunks per trans call: fewer holds less at once

    # discovery_run's pairs, in its order: the trans test fits each target with the same chunk-mates
    _passing = pairs_qc.loc[pairs_qc["pass_qc"], ["grna_target", "response_id"]]
    _reference = discovery_pairs[["response_id", "grna_target"]].merge(_passing, on=["grna_target", "response_id"])
    _cis = discovery_pairs.loc[discovery_pairs["in_fdr"].to_numpy(dtype=bool), ["response_id", "grna_target"]].merge(
        _passing, on=["grna_target", "response_id"]
    )

    # every element in the correction against every cis gene on another chromosome, chrY left out
    _targets = pd.unique(_cis["grna_target"])
    _genes = pd.unique(_cis["response_id"])
    _target_chr = np.repeat(elements.loc[_targets, "chr"].astype(str).to_numpy(), len(_genes))
    _gene_chr = np.tile(genes.loc[_genes, "gene_chr"].astype(str).to_numpy(), len(_targets))
    _keep = (_target_chr != _gene_chr) & (_gene_chr != "chrY")
    _trans_pairs = pd.DataFrame(
        {
            "response_id": np.tile(np.asarray(_genes), len(_targets))[_keep],
            "grna_target": np.repeat(np.asarray(_targets), len(_genes))[_keep],
        }
    )

    # discovery's target fit chunks, rebuilt from its pairs and its settings (run_discovery_analysis's defaults)
    _chunks = pysceptre.pipeline.trans.trans_fit_chunks(
        grna_target_cells, _reference, response_matrix.shape[1], target_chunk_size=200, chunk_memory_gb=1.0
    )
    _layout = (len(_chunks[0]), len(_chunks))  # (width, count): every call refuses any other layout
    _chunk_of = {_t: _i for _i, _chunk in enumerate(_chunks) for _t in _chunk}
    _trans_pairs["chunk"] = _trans_pairs["grna_target"].map(_chunk_of)
    _ids = np.sort(_trans_pairs["chunk"].unique())
    _calls = [_ids[_i : _i + TRANS_CHUNKS_PER_CALL] for _i in range(0, len(_ids), TRANS_CHUNKS_PER_CALL)]

    _parts = []
    for _call in mo.status.progress_bar(_calls, title="Trans pairs, a few fit chunks per call", show_eta=True):
        _parts.append(
            pysceptre.run_trans_analysis(
                response_matrix,
                gene_ids,
                covariate_matrix,
                grna_target_cells,
                _trans_pairs.loc[_trans_pairs["chunk"].isin(_call), ["response_id", "grna_target"]],
                reference_pairs=_reference,
                seed=SEED,
                side=SIDE,
                expected_fit_chunks=_layout,
                fit_chunk_ids=_call.tolist(),
                gene_fits=gene_fits,
                n_jobs=N_JOBS,
            )
        )
    trans_results = pd.concat(_parts, ignore_index=True)

    _screen = trans_results["screen"].value_counts()
    mo.md(
        f"**{len(trans_results):,}** trans pairs ({len(_targets):,} elements by {len(_genes)} genes) in "
        f"{len(_calls)} calls of up to {TRANS_CHUNKS_PER_CALL} of discovery's {_layout[1]} fit chunks "
        f"({_layout[0]} elements each): {_screen.get('none', 0):,} resampled, "
        f"{_screen.get('fold_change', 0):,} skipped for a fold change of at least 1, "
        f"{_screen.get('early_stop', 0):,} stopped once p was known to be at least 0.05."
    )
    return (trans_results,)


@app.cell
def long_range(
    FDR,
    SEED,
    discovery_table,
    elements,
    gene_ids,
    genes,
    mo,
    np,
    pd,
    pysceptre,
    trans_results,
):
    LONG_RANGE_BINS = [0, 10_000, 100_000, 500_000, np.inf]  # the class tables' bins, the long range split in two

    _in_fdr = discovery_table.loc[discovery_table["in_fdr"].fillna(False).to_numpy(dtype=bool)]
    _cutoff = pysceptre.analytical_power.bh_nominal_cutoff(_in_fdr["p_value"].to_numpy(), FDR)
    _elements = (
        elements.assign(centre=(elements["start"] + elements["end"]) / 2)
        .rename(columns={"chr": "chrom"})
        .rename_axis("grna_target")
        .reset_index()[["grna_target", "chrom", "centre"]]
    )
    _genes = genes.rename(columns={"gene_chr": "chrom", "gene_tss": "tss"}).rename_axis("response_id").reset_index()[
        ["response_id", "chrom", "tss"]
    ]

    long_range = pysceptre.run_specificity_check(
        _in_fdr[["grna_target", "response_id", "p_value", "fold_change"]],
        trans_results,
        _elements,
        _genes,
        _cutoff,
        measured_genes=list(gene_ids),
        distance_bins=LONG_RANGE_BINS,
        seed=SEED,
    )
    # the check's links are discovery's negative-effect hits, or the two disagree on the cutoff
    _negative = int(((_in_fdr["q_value"] < FDR) & (_in_fdr["fold_change"] < 1)).sum())
    assert int(long_range.links["link"].sum()) == _negative, (int(long_range.links["link"].sum()), _negative)


    def _ratios(by_distance):
        """The by-distance table as rates, with the trans/cis ratio (1 - above) and its interval."""
        _t = by_distance[["distance_bin", "pairs", "links"]].copy()
        _t["cis rate, %"] = (100 * by_distance["cis_rate"]).round(3)
        _t["trans level, %"] = (100 * by_distance["background"]).round(3)
        # a bin without links has no ratio: nothing there for the trans level to explain
        _some = (by_distance["links"] > 0).to_numpy()
        _t["trans/cis"] = (1 - by_distance["above"]).where(_some).round(2)
        _t["trans/cis, 95% interval"] = [
            f"{1 - _hi:.2f} to {1 - _lo:.2f}" if _s else ""
            for _s, _lo, _hi in zip(_some, by_distance["above_lo"], by_distance["above_hi"])
        ]
        return _t.set_index("distance_bin")


    _broad = int(long_range.broad_effect["broad_effect"].sum())
    _far = long_range.far_links["status"].value_counts() if "status" in long_range.far_links else pd.Series(dtype=int)
    mo.vstack(
        [
            mo.md(
                f"**{_negative}** links (negative effects at p ≤ {_cutoff:.2e}). The trans level: "
                f"**{long_range.background_rate:.3%}** of trans tests called. {_broad} elements have broad trans "
                "effects, far more calls than that level allows."
            ),
            _ratios(long_range.by_distance),
            mo.md(f"Without the {_broad} broad-effect elements, on both sides:"),
            _ratios(long_range.by_distance_without_broad),
            mo.md(
                "Long-range links (> 100 kb): could another gene the element lowers explain them, an indirect "
                "effect in cis? "
                + (", ".join(f"{_k} {_v}" for _k, _v in _far.items()) or "none")
                + "."
            ),
        ]
    )
    return LONG_RANGE_BINS, long_range


@app.cell
def long_range_plot(
    LONG_RANGE_BINS,
    links_vs_trans,
    long_range,
    mo,
    np,
    plot_links_vs_trans,
):
    _counts = links_vs_trans(long_range.by_distance)
    # the bins within 100 kb against those beyond it
    _upper = np.asarray(LONG_RANGE_BINS[1:], dtype=float)
    _near, _far = _counts[_upper <= 100_000], _counts[_upper > 100_000]


    def _share(d):
        _n = d["significant cis"].sum()
        return f"{d['cis above the trans level'].sum() / _n:.0%}" if _n else "none"


    _trans = long_range.background
    mo.vstack(
        [
            plot_links_vs_trans(
                _counts,
                title=f"{_share(_near)} of links within 100 kb clear the trans level, {_share(_far)} beyond",
                subtitle="Day 0: element-gene pairs significant with a negative effect, by distance to the TSS",
            ),
            mo.md(
                "Two bars per distance bin. Left, the significant cis pairs: in blue those above the trans "
                "level, in vermillion the number the trans level alone would make significant among the bin's "
                "tested pairs. Right, that trans-level number by itself: the bin's tested pairs times the "
                "trans rate of their genes. That number is an expected value, an average over repeats of "
                "the screen, so it can be a fraction: 0.5 means about even odds of one such pair. "
                "In the trans analysis itself "
                f"**{int(_trans['background'].sum()):,}** of {len(_trans):,} pairs are significant "
                f"({_trans['background'].mean():.4%})."
            ),
            _counts,
        ]
    )
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
def calibration_fns(np):
    def guide_cells(guide_type, guides, assigned):
        """Each guide of one type, with the positions of the cells it is assigned to."""
        calls = assigned.tocsr()
        rows = np.flatnonzero(guides["type"].astype(str).to_numpy() == guide_type)
        return {guides.index[i]: calls.indices[calls.indptr[i] : calls.indptr[i + 1]] for i in rows}

    return (guide_cells,)


@app.cell
def power_plot_fn(headline, house_style, np, plt):
    def plot_powered_by_effect(
        powered, totals=None, reference=None, reference_label="", legend_title="", title="", subtitle=""
    ):
        """Powered pairs against effect size, one line per ordered group, from none at no effect.

        `powered`: one row per effect size (0 to 1), one column per group, counts of pairs.
        `totals`: pairs per group; with `reference`, each group's dot there is labelled with its share.
        """
        # Okabe-Ito, warm to cool along the groups' order (yellow left out: too faint for a line on white),
        # and each one's text tint, at least 4.5:1 on white
        colours = ["#D55E00", "#E69F00", "#009E73", "#56B4E9", "#0072B2"]
        tints = ["#C35500", "#9F6D00", "#008661", "#0E7DB0", "#0072B2"]
        if powered.shape[1] > len(colours):
            raise ValueError(f"At most {len(colours)} groups, not {powered.shape[1]}.")
        effects = 100 * np.r_[0, powered.index.to_numpy(dtype=float)]
        top = powered.to_numpy().max() * 1.1
        fig, ax = plt.subplots(figsize=(7.5, 4.4))
        fig.subplots_adjust(left=0.11, right=0.95, bottom=0.12, top=0.74)
        if reference is not None:
            ax.axvline(100 * reference, color="0.1", linewidth=1)
            ax.annotate(reference_label, (100 * reference, 1), xycoords=("data", "axes fraction"),
                        xytext=(4, -9), textcoords="offset points", fontsize=8, color="0.1")
        # no pair is powered at no effect, where power is the test's own false-positive rate
        for (name, counts), colour in zip(powered.items(), colours):
            ax.plot(effects, np.r_[0, counts.to_numpy()] / 1e3, color=colour, linewidth=2, marker="o",
                    markersize=4, markevery=slice(1, None), label=name)
        if reference is not None and totals is not None:
            # below and right of its dot; dots on the axis stack up its left side instead, lowest group lowest
            stacked = 0
            for (name, counts), tint in zip(powered.items(), tints):
                count = counts.loc[reference]
                share = count / totals.loc[name]
                text = f"{share:.1%}" if 0 < share < 0.01 else f"{share:.0%}"
                if count < 0.04 * top:
                    offset, ha = (-5, 3 + 11 * stacked), "right"
                    stacked += 1
                else:
                    offset, ha = (5, -12), "left"
                ax.annotate(text, (100 * reference, count / 1e3), xytext=offset, textcoords="offset points",
                            ha=ha, fontsize=8.5, fontweight="bold", color=tint)
        ax.set_xticks(effects, [f"{e:g}" for e in effects])
        ax.set(xlim=(0, effects[-1] * 1.02), ylim=(0, top / 1e3), xlabel="Negative effect (%)", ylabel="Pairs, thousands")
        ax.tick_params(axis="x", length=0)
        ax.legend(loc="lower left", bbox_to_anchor=(0, 1.03), ncol=powered.shape[1], frameon=False,
                  fontsize=8.5, handlelength=1.6, columnspacing=1.4, title=legend_title, title_fontsize=8.5,
                  alignment="left", borderaxespad=0, borderpad=0)
        house_style(ax)
        headline(fig, title, subtitle)
        return fig


    return (plot_powered_by_effect,)


@app.cell
def hits_plot_fn(matplotlib, np, plt):
    def plot_hits_by_class(hits, tested=None, title="", subtitle=""):
        """Hits per class, negative effects leftward and positive rightward, each bar split by distance to the TSS.

        `hits`: one row per hit, with categorical `pair_class` and `distance` and boolean `down`.
        `tested`: pairs tested per class, printed under each class's name.
        """
        # Okabe-Ito, warm to cool from near to far, all at least 3:1 on white
        colours = ["#D55E00", "#009E73", "#0072B2"]
        gridline = "#B7C6CF"
        distances = list(hits["distance"].cat.categories)
        if len(distances) > len(colours):
            raise ValueError(f"At most {len(colours)} distance bins, not {len(distances)}.")
        counts = hits.groupby(["pair_class", "down", "distance"], observed=False).size()
        classes = list(hits.groupby("pair_class", observed=False).size().sort_values(kind="stable").index)  # most hits on top
        rows = np.arange(len(classes))
        height = 2.0 + 0.55 * len(classes)
        fig, ax = plt.subplots(figsize=(7.5, height))
        fig.subplots_adjust(left=0.2, right=0.95, bottom=0.5 / height, top=1 - 1.4 / height)
        widest = 0
        for side, down in ((-1, True), (1, False)):
            filled = np.zeros(len(classes))
            for distance, colour in zip(distances, colours):
                width = np.array([counts.loc[(c, down, distance)] for c in classes], dtype=float)
                ax.barh(rows, side * width, left=side * filled, height=0.62, color=colour, linewidth=0,
                        label=distance if down else None)
                filled += width
            widest = max(widest, filled.max())
            for row, total in zip(rows, filled):
                ax.annotate(f"{total:.0f}", (side * total, row), xytext=(4 * side, 0), textcoords="offset points",
                            ha="left" if side > 0 else "right", va="center", fontsize=8.5, color="0.1")
        limit = widest * 1.15
        ax.axvline(0, color="0.1", linewidth=1)
        ax.set_xlim(-limit, limit)
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{abs(v):g}"))
        names = [f"{c}\n{tested[c]:,} tested" if tested is not None else c for c in classes]
        ax.set_yticks(rows, names)
        ax.set_ylim(-0.6, len(classes) - 0.4)
        ax.set_xlabel("Hits")
        ax.set_axisbelow(True)
        ax.grid(axis="x", linewidth=0.5, color=gridline)
        ax.grid(axis="y", visible=False)
        ax.spines[["top", "right", "left", "bottom"]].set_visible(False)
        ax.tick_params(axis="both", length=0)
        ax.tick_params(axis="y", labelsize=8.5)
        for x, text, ha in ((-0.02, "← Negative effect", "right"), (0.02, "Positive effect →", "left")):
            ax.annotate(text, (0.5 + x, 1), xycoords="axes fraction", xytext=(0, 4), textcoords="offset points",
                        ha=ha, va="bottom", fontsize=9, fontweight="bold", color="0.1")
        ax.legend(loc="lower left", bbox_to_anchor=(-0.2 / 0.75 + 0.03 / 0.75, 1.1), ncol=len(distances), frameon=False,
                  fontsize=8.5, handlelength=1.2, columnspacing=1.4, title="Distance to the TSS", title_fontsize=8.5,
                  alignment="left", borderaxespad=0, borderpad=0)
        fig.text(0.03, 1 - 0.3 / height, title, fontsize=13, fontweight="bold", ha="left")
        fig.text(0.03, 1 - 0.55 / height, subtitle, fontsize=9.5, color="0.3", ha="left")
        return fig


    return (plot_hits_by_class,)


@app.cell
def long_range_plot_fn(headline, house_style, matplotlib, np, pd, plt):
    def links_vs_trans(by_distance):
        """Per distance bin, cis links and how many of them the trans level accounts for, as counts.

        `by_distance`: `run_specificity_check`'s `by_distance` table. The trans count of a bin is its tested
        pairs times its matched trans rate; the links above it are the rest, with the bootstrap interval.
        """
        expected = by_distance["pairs"] * by_distance["background"]
        above = (by_distance["links"] - expected).clip(lower=0)
        return pd.DataFrame(
            {
                "tested pairs": by_distance["pairs"].to_numpy(),
                "significant cis": by_distance["links"].to_numpy(),
                "significant trans, expected": expected.round(3).to_numpy(),
                "cis above the trans level": above.round(1).to_numpy(),
                "above, 95% interval": [
                    f"{max(n * lo, 0):.1f} to {n * hi:.1f}" if n > 0 else ""
                    for n, lo, hi in zip(by_distance["links"], by_distance["above_lo"].fillna(0), by_distance["above_hi"].fillna(0))
                ],
            },
            index=pd.Index(by_distance["distance_bin"].astype(str), name="distance to the TSS"),
        )


    def plot_links_vs_trans(counts, title="", subtitle=""):
        """Per distance bin, the cis links stacked as above and within the trans level, beside the trans count.

        `counts`: `links_vs_trans` output.
        """
        blue, vermillion, tint = "#0072B2", "#D55E00", "#C35500"  # Okabe-Ito; tint: vermillion text, 4.5:1
        x = np.arange(len(counts))
        links = counts["significant cis"].to_numpy(dtype=float)
        trans = counts["significant trans, expected"].to_numpy(dtype=float)
        within = np.minimum(trans, links)
        width = 0.34
        fig, ax = plt.subplots(figsize=(7.5, 4.4))
        fig.subplots_adjust(left=0.11, right=0.95, bottom=0.14, top=0.74)
        ax.bar(x - width / 2, within, width, color=vermillion, linewidth=0, label="Cis, within the trans level")
        ax.bar(x - width / 2, links - within, width, bottom=within, color=blue, linewidth=0, label="Cis, above the trans level")
        ax.bar(x + width / 2, trans, width, color=vermillion, linewidth=0, alpha=0.45,
               label="Trans level: significant pairs expected")
        top = max(links.max(), trans.max(), 1)
        for xi, n, t in zip(x, links, trans):
            ax.annotate(f"{n:.0f}", (xi - width / 2, n), xytext=(0, 3), textcoords="offset points", ha="center",
                        fontsize=9, fontweight="bold", color="0.1")
            label = "<0.1" if 0 < t < 0.05 else (f"{t:.1f}" if t < 10 else f"{t:.0f}")
            ax.annotate(label, (xi + width / 2, t), xytext=(0, 3),
                        textcoords="offset points", ha="center", fontsize=8.5, color=tint)
        ax.set_xticks(x, [f"{b}\n{n:,} tested" for b, n in zip(counts.index, counts["tested pairs"])])
        ax.set(xlim=(-0.6, len(counts) - 0.4), ylim=(0, top * 1.15), ylabel="Significant pairs")
        ax.tick_params(axis="x", length=0, labelsize=8.5)
        ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
        handles, labels = ax.get_legend_handles_labels()
        order = [1, 0, 2]  # blue first: the cis links above the trans level are the point
        ax.legend([handles[i] for i in order], [labels[i] for i in order], loc="lower left", bbox_to_anchor=(0, 1.03), ncol=3, frameon=False, fontsize=8.5, handlelength=1.2,
                  columnspacing=1.4, borderaxespad=0, borderpad=0)
        house_style(ax)
        headline(fig, title, subtitle)
        return fig


    return links_vs_trans, plot_links_vs_trans


if __name__ == "__main__":
    app.run()
