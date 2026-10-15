# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "altair>=6.0.0",
#     "anndata==0.13.4",
#     "crispr-analysis-utils",
#     "h5py>=3.12.0",
#     "marimo>=0.20.2",
#     "matplotlib==3.11.2",
#     "mudata==0.4.1",
#     "numpy>=2.1.0",
#     "pandas>=2.2.0",
#     "scipy>=1.14.0",
# ]
#
# [tool.uv.sources]
# crispr-analysis-utils = { git = "https://github.com/EngreitzLab/crispr-analysis-utils", rev = "6809bb5521220bf2c74c758e131f9e4398436ba3" }
# ///

import marimo

__generated_with = "0.25.1"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell
def _():
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp

    from crispr_analysis_utils import qc


    return mpl, np, pd, plt, qc


@app.cell
def _(mpl, qc):
    # The house style: Okabe-Ito on the colour cycle, cividis on images, a grey
    # grid behind the data, no top or right spine. Held in the package so a
    # figure drawn here and one drawn by the pipeline look the same.
    qc.style()

    # Screen only: the published PDF is written at qc.DPI regardless.
    mpl.rcParams["figure.dpi"] = 140
    mpl.rcParams["savefig.facecolor"] = "white"

    return


@app.cell
def _(mo):
    from pathlib import Path

    # Everything resolves from the repository, so no path here names a machine.
    REPO = mo.notebook_dir().parents[1]

    # Data is gitignored in full, so nothing large is ever tracked.
    DATA = REPO / "data"

    # The pipeline's per-channel kb output, one folder per measurement set.
    PER_CHANNEL = "gs://hgrm-seqera/igvf-crispr-all/per-channel"

    # Where a channel is pulled down to before it is read.
    LOCAL = DATA / "per_channel"
    LOCAL.mkdir(parents=True, exist_ok=True)

    # Figures, their numbers and their captions, beside this notebook.
    RESULTS = mo.notebook_dir() / "results"

    # The lab's Cell Ranger output for the channel the pipeline was tested on.
    CELLRANGER = (
        DATA
        / "h5"
        / "rep01_day04_chip_a_batch01_channel03_sample_filtered_feature_bc_matrix.h5"
    )

    return DATA, Path, REPO, RESULTS


@app.cell
def _(mo):
    mo.md("""
    # Pipeline counts against Cell Ranger

    Per-channel kb output from the IGVF CRISPR pipeline, read straight from
    the staged copy rather than through a merged MuData.
    """)
    return


@app.cell
def _(DATA):
    # The assembled object: one channel, both modalities, nothing filtered.
    MUDATA = DATA / "mudata" / "mudata.h5mu"
    # What this notebook writes back: the same object, QC-passing cells only.
    QC_MUDATA = DATA / "mudata" / "mudata_qc.h5mu"

    return MUDATA, QC_MUDATA


@app.cell
def _(MUDATA):
    import mudata as md
    from mudata._core.settings import settings as _mudata_settings

    # mudata defaults display_style to "text", and in that mode _repr_html_
    # returns None, so a cell ending in the object renders nothing at all.
    _mudata_settings.display_style = "html"

    mdata = md.read_h5mu(MUDATA)
    mdata
    return (mdata,)


@app.cell
def _(mdata):
    mdata["gene"].obs
    return


@app.cell
def _(mdata):
    # mudata's HTML repr summarises the slots and never lists var columns,
    # whatever display_html_expand is set to, so look at the tables themselves.
    mdata["gene"].var
    return


@app.cell
def _(mdata):
    mdata["guide"].var
    return


@app.cell
def _(mdata):
    # One row per guide per site per gene, so it lives in uns rather than var.
    mdata["guide"].uns["annotation"]
    return


@app.cell
def _(mdata):
    mdata["guide"].uns["element_mapping"]
    return


@app.cell
def _(mdata, mo, qc):
    # A cell missing from a modality has no UMIs there either, so both count.
    gene_umis = qc.umis_per_cell(mdata, "gene")
    guide_umis = qc.umis_per_cell(mdata, "guide")

    no_gene_umis = int((gene_umis == 0).sum())
    no_guide_umis = int((guide_umis == 0).sum())

    mo.md(f"""
    Of {mdata.n_obs:,} cells:

    - **{no_gene_umis:,}** have no gene UMIs
    - **{no_guide_umis:,}** have no guide UMIs
    """)

    return


@app.cell
def _(RESULTS, mdata, qc):
    # The object names its own channel, so nothing here is hardcoded.
    CHANNEL = f"{mdata.uns['channel']['channel_name']} ({mdata.uns['measurement_set']})"

    # qc.knee ranks every barcode in both modalities and returns the figure, the
    # numbers behind it and a caption generated from the same counts.
    knee_figure, knee_numbers, knee_caption = qc.knee(mdata)
    knee_figure.axes[0].set_title(f"Knee plot, {CHANNEL}")

    qc.publish(RESULTS / "knee_umis_per_barcode", knee_figure, knee_numbers, knee_caption)
    knee_figure

    return (CHANNEL,)


@app.cell
def _(mdata, mo, qc):
    # A first cut, read off the knee above: below this a barcode is ambient.
    # Applied to the gene UMIs, so a barcode kept here is one with real RNA.
    UMI_FLOOR = 1000

    cells, keep = qc.filter_cells(mdata, min_umis=UMI_FLOOR, modality="gene")

    mo.md(f"""
    At **{UMI_FLOOR:,}** gene UMIs: kept **{int(keep.sum()):,}** of {mdata.n_obs:,} barcodes.
    """)

    return UMI_FLOOR, cells


@app.cell
def _(cells, qc):
    # scanpy names these as the IGVF pipeline does, so they are directly
    # comparable. mito is the chrM contig, ribo is RPS/RPL on protein-coding
    # genes, ribo_pseudogene the same prefixes elsewhere, targeted the TAP panel.
    #
    # ribo_rmd is this screen's own convention, from the R QC
    # (analyses/2026-01-20_load_counts_and_qc): wider than the module's, and the
    # definition the per-day ribosomal thresholds were read against. Seurat's
    # PercentageFeatureSet greps the symbols, so this does too.
    RIBO_RMD = r"^RPS|^RPL|^RRP|^MRP|^FAU"
    ribo_rmd = (
        cells["gene"].var["symbol"].astype("string").fillna("")
        .str.contains(RIBO_RMD, regex=True)
        .to_numpy()
    )

    fractions = qc.cell_fractions(cells, modality="gene", extra={"ribo_rmd": ribo_rmd})
    fractions.describe().T[["mean", "50%", "min", "max"]].round(2)

    return (fractions,)


@app.cell
def _(mdata, mo):
    # Specific to this screen: the thresholds Eugenio settled on per
    # differentiation day in the R QC (analyses/2026-01-20_load_counts_and_qc),
    # where they were read off the day's own distributions rather than derived.
    # min_ribo is on the R definition of ribosomal, which is wider than ours.
    DAY_THRESHOLDS = {
        "0": {"min_umis": 1500, "min_targeted": 55, "min_ribo": 20},
        "2": {"min_umis": 1000, "min_targeted": 45, "min_ribo": 16},
        "4": {"min_umis": 1000, "min_targeted": 25, "min_ribo": 6},
    }

    DAY = mdata.uns["channel"]["day"]
    CUTS = DAY_THRESHOLDS[DAY]

    mo.md(f"Day **{DAY}**: {CUTS}")

    return CUTS, DAY


@app.cell
def _(CUTS, fractions, mo, qc):
    # Two kinds of rule. A band measured from this channel's own spread, for the
    # three metrics where a cell far from the middle is a cell gone wrong; and
    # this screen's day thresholds, read off the day's own distributions in the R
    # QC rather than derived. Nothing is applied yet -- this is the list.
    NMADS = 5
    BAND = ("total_counts", "n_genes_by_counts", "pct_counts_mito")
    LIMITS = {
        "total_counts": (CUTS["min_umis"], None),
        "pct_counts_targeted": (CUTS["min_targeted"], None),
        "pct_counts_ribo_rmd": (CUTS["min_ribo"], None),
    }

    passes, rules = qc.cell_filters(fractions, nmads=NMADS, band=BAND, limits=LIMITS)

    mo.md(f"""
    Of {len(passes):,} cells, each rule on its own:

    {rules.to_markdown(index=False, floatfmt=",.1f")}
    """)

    return BAND, LIMITS, NMADS, passes, rules


@app.cell
def _(BAND, CHANNEL, LIMITS, NMADS, RESULTS, fractions, qc):
    # The same rules the filter applies, so the lines here are the cuts the
    # written object was made with: a band on the three metrics in BAND, the
    # screen's day thresholds on the rest. Blue passes every rule a metric
    # carries, grey fails one.
    dist_figure, dist_numbers, dist_caption = qc.metric_distributions(
        fractions, nmads=NMADS, band=BAND, limits=LIMITS
    )
    dist_figure.suptitle(f"QC metrics, {CHANNEL}", y=1.0)

    qc.publish(RESULTS / "qc_metric_distributions", dist_figure, dist_numbers, dist_caption)
    dist_figure

    return


@app.cell
def _(CHANNEL, CUTS, RESULTS, fractions, qc):
    # The R analysis's central QC plot: depth against how much of a cell is
    # the readout panel. A cell can fail on either alone.
    scatter_figure, scatter_numbers, scatter_caption = qc.metric_scatter(
        fractions,
        x="total_counts",
        y="pct_counts_targeted",
        x_threshold=CUTS["min_umis"],
        y_threshold=CUTS["min_targeted"],
    )
    scatter_figure.suptitle(f"Depth against readout share, {CHANNEL}")

    qc.publish(
        RESULTS / "total_umi_vs_targeted_percent",
        scatter_figure,
        scatter_numbers,
        scatter_caption,
    )
    scatter_figure

    return


@app.cell
def _(QC_MUDATA, REPO, UMI_FLOOR, cells, mdata, mo, passes, qc, rules):
    # Every rule at once, and the object written with the thresholds that made
    # it. A cell fails each rule independently, so the fails do not add up to
    # the cells lost: one cell can fail several.
    final = qc.apply_filters(cells, passes, rules, path=QC_MUDATA)

    _steps = "\n".join(
        f"| {_r.metric} {_r.source} | {_r.n_failed:,} | |" for _r in rules.itertuples()
    )

    mo.md(f"""
    | step | fails | left |
    | --- | ---: | ---: |
    | barcodes | | {mdata.n_obs:,} |
    | gene UMIs >= {UMI_FLOOR:,} (knee) | {mdata.n_obs - cells.n_obs:,} | {cells.n_obs:,} |
    {_steps}
    | **every rule** | | **{final.n_obs:,}** |

    Written to `{QC_MUDATA.relative_to(REPO)}`, every threshold above recorded in
    its `uns["qc"]`.
    """)

    return (final,)


@app.cell
def _(DAY, RESULTS, final, qc):
    # Does the panel dominate, and across the whole panel or on a few genes?
    panel_figure, panel_numbers, panel_caption = qc.panel_expression(final)
    panel_figure.axes[0].set_title(f"Panel vs background, day {DAY}", fontsize="medium")

    qc.publish(
        RESULTS / "targeted_vs_background_upm", panel_figure, panel_numbers, panel_caption
    )
    panel_figure

    return (panel_numbers,)


@app.cell
def _(CHANNEL, DATA, DAY, RESULTS, final, np, panel_numbers, pd, plt, qc):
    # Specific to this screen: the panel carries marker genes for each
    # differentiation day plus housekeeping genes, so their relative expression
    # says whether this channel sits where its day label claims.
    PANEL_META = DATA.parent / "metadata" / "metadata_genes_final.csv"

    def unversioned(ids):
        """GENCODE's own rule: drop a trailing .N, never a PAR_Y suffix."""
        return pd.Series(ids, dtype="string").str.replace(r"\.\d+$", "", regex=True)

    markers = (
        pd.read_csv(PANEL_META)
        .assign(gene=lambda d: unversioned(d["intended_target_gene_id"]))
        .merge(
            panel_numbers.assign(gene=lambda d: unversioned(d["gene"]))[["gene", "mean_upm"]],
            on="gene",
            how="left",
        )
    )
    # One row per gene: four genes carry two primer pairs and would count twice.
    markers = markers.drop_duplicates("gene")

    CLASSES = {
        "is_d0_marker": "Day 0",
        "is_d2_marker": "Day 2",
        "is_d4_marker": "Day 4",
        "is_housekeeping_gene": "Housekeeping",
    }
    marker_rows = pd.concat(
        [
            markers.loc[markers[column].astype(bool), ["gene", "intended_target_gene_symbol", "mean_upm"]].assign(gene_class=label)
            for column, label in CLASSES.items()
        ],
        ignore_index=True,
    )

    # One label per panel gene, in the R's own precedence order, so the same
    # classes can colour any per-gene figure. A gene in no class stays NA.
    def _class_of(row):
        for _column, _label in CLASSES.items():
            if bool(row[_column]):
                return _label
        return None

    panel_class = pd.Categorical(
        markers.set_index("gene").apply(_class_of, axis=1),
        categories=list(CLASSES.values()),
    )
    panel_class = pd.Series(panel_class, index=markers["gene"])

    _fig, _ax = plt.subplots(figsize=(5.0, 4.0))
    _rng = np.random.default_rng(0)
    for _i, _label in enumerate(CLASSES.values()):
        _v = marker_rows.loc[marker_rows["gene_class"] == _label, "mean_upm"].to_numpy()
        # Outlined, like the per-gene scatter: the yellow is a pale fill.
        _ax.scatter(_i + _rng.uniform(-0.18, 0.18, len(_v)), _v, s=14,
                    color=qc.OKABE_ITO[_i], linewidths=qc.EDGE_WIDTH,
                    edgecolors=qc.EDGE)
        _ax.hlines(np.median(_v), _i - 0.3, _i + 0.3, color="black", linewidth=1.5)
    _ax.set(xticks=range(len(CLASSES)), xticklabels=list(CLASSES.values()),
            yscale="log", ylabel="mean UMIs per million",
            title=f"Marker expression, day {DAY}")
    _ax.grid(axis="y", which="major", color="0.85", linewidth=0.6)
    _ax.set_axisbelow(True)
    _fig.tight_layout()

    marker_caption = (
        f"Mean expression of the panel's marker and housekeeping genes in "
        f"{final.n_obs:,} cells of {CHANNEL}, one point per gene, black bars at "
        "the class medians. The day-4 markers should stand above the day-0 ones "
        "in a day-4 channel; the housekeeping genes are the scale to read them "
        "against. Medians: "
        + ", ".join(
            f"{_l} {marker_rows.loc[marker_rows['gene_class'] == _l, 'mean_upm'].median():,.0f}"
            for _l in CLASSES.values()
        )
        + "."
    )
    qc.publish(RESULTS / "marker_genes_upm", _fig, marker_rows, marker_caption)
    _fig

    return panel_class, unversioned


@app.cell
def _(DAY, Path, final, mo, np, pd, qc, unversioned):
    # Specific to this screen: the 10x multiome time course is the untargeted
    # reference. A panel gene well expressed there but unseen here is a primer
    # that did not work; one seen here but not there is the panel reaching
    # something the whole-transcriptome assay misses.
    REFERENCE_RNA = (
        Path.home()
        / "GitHub/broad-nnfc-cmo-multiome/results/10X_multiome_5_timepoints/merged/rna"
        / f"rna_d{DAY}.h5ad"
    )

    import anndata as ad

    _ref = ad.read_h5ad(REFERENCE_RNA)
    # Already strict-QC'd singlets; keep the CMO-tagged ones, since the
    # Negatives carry no timepoint and would mix days into a day reference.
    _ref = _ref[_ref.obs["timepoint_scanpy"] == f"d{DAY}"].copy()
    _ref.X = _ref.layers["counts"]
    reference_upm = pd.DataFrame(
        {"gene": _ref.var["gene_id_base"].to_numpy(), "reference_upm": qc.mean_upm(_ref)}
    ).groupby("gene", as_index=False)["reference_upm"].max()

    # Ours: mean UMIs per gene per cell, unnormalised, as the R analysis had it.
    _ours = final["gene"]
    _panel = _ours.var["readout"].fillna(False).to_numpy().astype(bool)
    ours_umi = pd.DataFrame(
        {
            "gene": unversioned(_ours.var_names[_panel]),
            "symbol": _ours.var["symbol"].to_numpy()[_panel],
            "mean_umi": np.asarray(_ours.X[:, _panel].mean(axis=0)).ravel(),
        }
    )

    comparison = ours_umi.merge(reference_upm, on="gene", how="left").set_index("gene")
    mo.md(
        f"**{len(comparison):,}** panel genes; "
        f"**{int(comparison['reference_upm'].isna().sum())}** absent from the "
        f"{_ref.n_obs:,}-cell day-{DAY} reference."
    )

    return (comparison,)


@app.cell
def _(DAY, RESULTS, comparison, panel_class, qc):
    # scRNA against TAP, one point per panel gene: how much of a gene the
    # untargeted multiome sees against how much this targeted assay sees.
    # Nothing is filtered here, so colour carries the gene's class and the
    # lines only mark where the R analysis expected the panel to fall. Those
    # thresholds, 30 UPM and 1 mean UMI per cell, were tuned on the 5' pilot.
    REF_UPM_CUT, OURS_UMI_CUT = 30, 1

    ref_figure, ref_numbers, ref_caption = qc.metric_scatter(
        comparison,
        x="reference_upm",
        y="mean_umi",
        x_threshold=REF_UPM_CUT,
        y_threshold=OURS_UMI_CUT,
        log=("reference_upm", "mean_umi"),
        unit="genes",
        groups=panel_class.reindex(comparison.index),
    )
    ref_figure.axes[0].set_xlabel("scRNA: mean UPM in the multiome reference")
    ref_figure.axes[0].set_ylabel("TAP: mean UMIs per cell")
    ref_figure.suptitle(f"scRNA vs TAP, {len(comparison)} panel genes, day {DAY}")

    # Carry the symbol through, so a gene in a corner can be named.
    ref_numbers = ref_numbers.merge(
        comparison.reset_index()[["gene", "symbol"]], on="gene", how="left"
    )
    qc.publish(RESULTS / "scrna_vs_tap", ref_figure, ref_numbers, ref_caption)
    ref_figure

    return


@app.cell
def _(final, mo, qc):
    # Guide QC, on the cells that passed. The guide modality holds fewer cells
    # than the object does, so every guide's UMIs and the cells it reached are
    # summed over its own.
    guide_totals = qc.guide_counts(final)

    mo.md(f"""
    {guide_totals.attrs["n_cells"]:,} of the {final.n_obs:,} cells that passed QC
    carry the guide modality. They hold {int(guide_totals["umis"].sum()):,} guide
    UMIs over {len(guide_totals):,} guides, of which
    {int((guide_totals["umis"] == 0).sum()):,} appear in none of them.
    """)
    return (guide_totals,)


@app.cell
def _(CHANNEL, RESULTS, guide_totals, qc):
    # How much of the library this channel saw. The grey bar is the guides no
    # cell carried, which the log axis has no place for; ten thousand cells
    # cannot sample a library of forty thousand, so most of it is absent by
    # arithmetic rather than by a defect.
    guide_figure, guide_numbers, guide_caption = qc.guide_abundance(guide_totals)
    guide_figure.suptitle(f"Guide abundance, {CHANNEL}", y=1.0)

    qc.publish(RESULTS / "guide_abundance", guide_figure, guide_numbers, guide_caption)
    guide_figure
    return


if __name__ == "__main__":
    app.run()
