# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "anndata==0.13.4",
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
    # Day-0 discovery with pysceptre

    Day 0 of the hiPSC-to-endothelial-cell differentiation screen, read out on a
    targeted panel of 292 genes: which elements regulate which genes. It starts
    from the cells the day-0 QC notebook kept (`analyses/2026-11-23_day_objects_qc`)
    and rebuilds the guide calls, the pairs and the covariates exactly as that
    notebook did. Then come the pysceptre tutorial's calibration check, power
    check, discovery, power and long-range steps
    (`analyses/2026-11-17_pysceptre_tutorial`).

    The six long cells start disabled: `calibration_nt`, `calibration_safe`,
    `power_check`, `discovery_run`, `gene_fits` and `trans_results`. Enable them
    in that order.
    """)
    return


@app.cell
def import_marimo():
    import marimo as mo

    return (mo,)


@app.cell(hide_code=True)
def imports(mo):
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
    return matplotlib, mudata, np, pd, plt, pysceptre, sparse, stats


@app.cell(hide_code=True)
def md_get_data(mo):
    mo.md(r"""
    ## Get the data

    The cells `day0_qc.py` kept and wrote: `pass_qc` and the R QC's three day-0
    limits, with `gene_tss` and `gene_strand` on the genes and the limits in
    `uns["cell_limits"]`.
    """)
    return


@app.cell
def get_data(mo):
    data_path = mo.notebook_dir().parents[1] / "results" / "2026-11-23_day_objects_qc" / "day0_pass_qc.h5mu"
    mo.md(f"`{data_path}`: {data_path.stat().st_size / 1e9:.2f} GB")
    return (data_path,)


@app.cell(hide_code=True)
def md_load(mo):
    mo.md(r"""
    ## Load the screen

    Two modalities over the same cells: `gene` (UMIs per gene) and `guide` (UMIs
    per guide, plus a `guide_assignment` layer saying which guides each cell
    carries).
    """)
    return


@app.cell
def load(data_path, mo, mudata):
    mdata = mudata.read_h5mu(data_path)
    mo.vstack(
        [
            mo.md(f"**{mdata.n_obs:,}** cells; limits {dict(mdata.uns['cell_limits'])}"),
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
def guide_threshold(mo):
    guide_threshold = mo.ui.slider(
        start=20, stop=200, step=10, value=20, label="Guide UMI threshold", show_value=True
    )
    guide_threshold
    return (guide_threshold,)


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
    sides: at least 20 **treated** cells (carrying any of the element's guides) and
    90,000 **control** cells (all the other cells) with a nonzero count. sceptre's
    defaults are 7 and 7. The control cut is raised to drop the genes day 0 barely
    expresses: the readout genes split at about 100,000 nonzero cells, as the
    scatter below shows. pysceptre's power check applies this itself; discovery
    expects pairs that already pass, so we apply it here, with the same function,
    to every pair.
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
                "all pass: day 0 barely expresses them."
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
                subtitle="Day 0: cells with a nonzero count of the pair's gene; log scales",
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
    return covariate_matrix, covariate_names


@app.cell(hide_code=True)
def md_gene_fit_check(mo):
    mo.md(r"""
    ### Check: are the gene fits degenerate?

    The first calibration run warned that 80 of 208 gene fits had a numerically
    rank-deficient design: the smallest eigenvalue of the weighted covariate
    cross-product sits within roundoff of zero, so statistics built from it are
    unreliable. pysceptre flags a fit when that eigenvalue is at most the largest
    × the number of covariates × machine epsilon.

    Twelve tested genes, spread across the expression range, fitted four ways:
    the current design (five covariates and 35 channel dummies), without the
    channel dummies, with replicate dummies instead, all in one job; and the
    current design again with `N_JOBS` workers, as the calibration fits them.

    On this Mac pysceptre runs parallel jobs as threads (Linux forks processes).
    The threaded run gives the same coefficients and theta, bit for bit, but
    negative eigenvalues, which a weighted cross-product cannot have: the
    precomputation is not safe to run in several threads at once. The degenerate
    fits are that, not the design.
    """)
    return


@app.cell
def gene_fit_check(
    N_JOBS,
    cell_ids,
    covariate_matrix,
    covariate_names,
    gene_ids,
    genes,
    mdata,
    mo,
    np,
    pairs_qc,
    pd,
    pysceptre,
    response_matrix,
):
    # twelve genes of the pairs that pass QC, evenly spaced by mean UMIs per cell
    _mean = pd.Series(np.asarray(response_matrix.sum(axis=1)).ravel() / response_matrix.shape[1], index=gene_ids)
    _tested = _mean[pd.unique(pairs_qc.loc[pairs_qc["pass_qc"], "response_id"])].sort_values()
    _pick = _tested.index[np.linspace(0, len(_tested) - 1, 12).round().astype(int)]
    _rows = pd.Index(gene_ids).get_indexer(_pick).tolist()

    _base = [_i for _i, _n in enumerate(covariate_names) if not _n.startswith("channel_")]
    _replicate = pd.get_dummies(mdata.obs.loc[cell_ids, "replicate"].astype(str), drop_first=True).to_numpy(dtype=float)
    _designs = {
        "channel": covariate_matrix,
        "none": covariate_matrix[:, _base],
        "replicate": np.hstack([covariate_matrix[:, _base], _replicate]),
    }

    # the calibration's own call: one shared outer-product array, N_JOBS workers
    _runs = {_name: (_z, 1) for _name, _z in _designs.items()} | {f"channel, {N_JOBS} jobs": (covariate_matrix, N_JOBS)}

    _rows_out = []
    for _name, (_z, _jobs) in _runs.items():
        _fits = pysceptre.pipeline.discovery.fit_all_genes(
            response_matrix,
            list(_pick),
            _z,
            gene_rows=_rows,
            batch_width=1,
            n_jobs=_jobs,
            x_outer_flat_shared=pysceptre.glm.irls.x_outer_flat(_z) if _jobs > 1 else None,
        )
        _flags = pysceptre.pipeline.discovery.summarize_gene_fits(_fits)
        for _g, _fit in _fits.items():
            _rows_out.append(
                {
                    "design": _name,
                    "covariates": _z.shape[1],
                    "jobs": _jobs,
                    "gene": genes.loc[_g, "symbol"],
                    "mean_umi": _mean[_g],
                    "eig_ratio": _fit.min_eigenvalue / _fit.max_eigenvalue,
                    "singular": _g in _flags["singular_design"],
                    "glm_converged": bool(_fit.glm_converged),
                    "theta_fallback": _g in _flags["theta_fallback"],
                }
            )
    gene_fit_check = pd.DataFrame(_rows_out)

    _singular = gene_fit_check.groupby("design", sort=False)["singular"].sum()
    mo.vstack(
        [
            mo.md(
                "Singular fits of 12: "
                + ", ".join(f"**{_d}** {int(_n)}" for _d, _n in _singular.items())
                + ". `eig_ratio` is the smallest eigenvalue over the largest."
            ),
            gene_fit_check.pivot(index=["gene", "mean_umi"], columns="design", values=["eig_ratio", "singular"]).round(
                {"mean_umi": 2}
            ),
        ]
    )
    return


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
    N_JOBS = 8  # parallel jobs

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

    _rows, _per_hit = [], []
    for _pair in _hits.itertuples():
        _gene = readout_genes.loc[_pair.response_id, "symbol"]
        _expression = response_matrix[gene_ids.index(_pair.response_id)].toarray().ravel() / _depth
        _guides = _pair.grna_target.split("&")
        # a cell can carry several of the element's guides: count it once
        _element_cells = len(np.unique(np.concatenate([control_guide_cells["non-targeting"][_g] for _g in _guides])))
        _fold = {_g: _expression[control_guide_cells["non-targeting"][_g]].mean() / _expression.mean() for _g in _guides}
        # the guide furthest from 1 on a log scale
        _top = max(_guides, key=lambda _g: abs(np.log2(max(_fold[_g], 1e-3))))
        _hit = {"gene": _gene, "pair_fold_change": _pair.fold_change, "p_value": _pair.p_value, "element_cells": _element_cells}
        _per_hit.append({**_hit, "guides": _guides, "top_guide": _top, "top_guide_fold_change": _fold[_top]})
        _rows += [{**_hit, "guide": _g, "guide_fold_change": _fold[_g]} for _g in _guides]
    calibration_hit_guides = pd.DataFrame(
        _rows, columns=["gene", "pair_fold_change", "p_value", "element_cells", "guide", "guide_fold_change"]
    )

    mo.vstack(
        [
            mo.md(f"**{len(_hits)}** non-targeting pairs pass Benjamini-Hochberg at {FDR}. Per hit, the element's guides and the one that moves the gene most:"),
            pd.DataFrame(_per_hit),
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
    # p < 0.05 uncorrected, as in sceptre's power check; empty where the pair failed QC and was not tested
    power_check.insert(
        power_check.columns.get_loc("p_value") + 1,
        "significant",
        power_check["p_value"].lt(0.05).astype("boolean").where(power_check["pass_qc"]),
    )
    _tested = power_check[power_check["pass_qc"]]
    _knockdown = (_tested["significant"] & (_tested["fold_change"] < 1)).astype(bool)

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


@app.cell(hide_code=True)
def md_discovery_ihw(mo):
    mo.md(r"""
    ### IHW, weighted by distance to the TSS

    Independent hypothesis weighting learns a weight for each stratum of a covariate
    from the data, then applies a weighted Benjamini-Hochberg at the same FDR. The
    covariate must be independent of the p-value under the null. Power is not: the
    safe-targeting calibration's p-values run small for highly expressed genes.
    Distance to the TSS should be, since the test never sees where an element is.

    IHW runs twice: with the 20 strata pyIHW chooses for these pairs, kept for
    reference, and with 5, whose weights agree better across folds (see "Fewer
    strata" below).
    """)
    return


@app.cell
def discovery_ihw(FDR, SEED, discovery_table, mo, np):
    import pyihw

    IHW_STRATA = (20, 5)  # 20 is pyIHW's own choice for these pairs; 5 steadies the weights (ihw_bins)

    # the pairs in the correction, with the distance to the TSS as the covariate
    _pairs = discovery_table.loc[
        discovery_table["in_fdr"].fillna(False).to_numpy(dtype=bool),
        ["grna_target", "response_id", "gene_symbol", "genomic_element", "distance_to_tss", "fold_change", "p_value", "q_value"],
    ].reset_index(drop=True)
    _bh = _pairs["q_value"] < FDR

    ihw_fits, ihw_strata, _columns, _hits = {}, {}, {}, {}
    for _n in IHW_STRATA:
        _fit = pyihw.ihw(
            _pairs["p_value"].to_numpy(), _pairs["distance_to_tss"].abs().to_numpy(), FDR, nbins=_n, seed=SEED
        )
        _hits[_n] = _fit.adj_pvalues() < FDR
        _columns |= {
            f"ihw{_n}_stratum": _fit.groups_factor(),
            f"ihw{_n}_weight": _fit.weights(),
            f"ihw{_n}_q_value": _fit.adj_pvalues(),
        }
        ihw_strata[_n] = (
            _pairs.assign(
                distance_kb=_pairs["distance_to_tss"].abs() / 1000,
                stratum=_fit.groups_factor(),
                weight=_fit.weights(),
                bh_hit=_bh,
                ihw_hit=_hits[_n],
            )
            .groupby("stratum")
            .agg(
                from_kb=("distance_kb", "min"),
                to_kb=("distance_kb", "max"),
                pairs=("p_value", "size"),
                mean_weight=("weight", "mean"),
                bh_hits=("bh_hit", "sum"),
                ihw_hits=("ihw_hit", "sum"),
            )
            .round(2)
        )
        # its rows are the rows of the fit's weights(levels_only=True)
        assert list(ihw_strata[_n].index) == list(range(1, _n + 1))
        ihw_fits[_n] = _fit
    discovery_ihw = _pairs.assign(**_columns)

    _any = np.logical_or.reduce(list(_hits.values()))
    mo.vstack(
        [
            mo.md(
                f"**{int(_bh.sum())}** hits at FDR {FDR} with Benjamini-Hochberg. With IHW: "
                + "; ".join(
                    f"**{int(_h.sum())}** with {_n} strata ({int((_h & ~_bh).sum())} gained, {int((_bh & ~_h).sum())} lost)"
                    for _n, _h in _hits.items()
                )
                + ". A weight above 1 lowers the bar for its stratum's pairs."
            ),
            *[mo.vstack([mo.md(f"**{_n} strata**"), ihw_strata[_n]]) for _n in IHW_STRATA],
            mo.ui.table(
                discovery_ihw.loc[_bh | _any].sort_values(f"ihw{IHW_STRATA[0]}_q_value"), selection=None
            ),
        ]
    )
    return IHW_STRATA, discovery_ihw, ihw_fits, ihw_strata, pyihw


@app.cell(hide_code=True)
def md_ihw_qc(mo):
    mo.md(r"""
    #### Checking the weighting

    IHW learns the weights in five folds, each from the pairs outside it, so folds that
    disagree mean a weight the data barely supports. Above a p-value of 0.5 nearly every
    pair is null, so each stratum's p-values should be flat there: a tilt would mean
    distance changes the null, and the FDR would no longer be controlled. Last, the hits
    per stratum under each correction. Each check has two plots: 20 strata, then 5.
    """)
    return


@app.cell
def ihw20_weights_plot(ihw_fits, ihw_strata, plot_ihw_weights):
    _n = 20
    _labels = [f"{_kb:.0f}" for _kb in ihw_strata[_n]["to_kb"]]
    _weights = ihw_fits[_n].weights(levels_only=True)  # strata x folds, the strata in ihw_strata's order
    plot_ihw_weights(
        _weights,
        _labels,
        xlabel="Distance to the TSS, kb: each stratum's upper edge",
        title=f"Pairs within {_labels[0]} kb get {_weights[0].mean():.1f} times the weight",
        subtitle=f"Day 0: IHW's weight per stratum of distance to the TSS, {_n} strata; Benjamini-Hochberg's is 1",
    )
    return


@app.cell
def ihw5_weights_plot(ihw_fits, ihw_strata, plot_ihw_weights):
    _n = 5
    _labels = [f"{_kb:.0f}" for _kb in ihw_strata[_n]["to_kb"]]
    _weights = ihw_fits[_n].weights(levels_only=True)  # strata x folds, the strata in ihw_strata's order
    plot_ihw_weights(
        _weights,
        _labels,
        xlabel="Distance to the TSS, kb: each stratum's upper edge",
        title=f"Pairs within {_labels[0]} kb get {_weights[0].mean():.1f} times the weight",
        subtitle=f"Day 0: IHW's weight per stratum of distance to the TSS, {_n} strata; Benjamini-Hochberg's is 1",
    )
    return


@app.cell
def ihw20_pvalues_plot(
    discovery_ihw,
    ihw_strata,
    np,
    plot_pvalues_by_stratum,
    stats,
):
    _n = 20
    _strata = discovery_ihw[f"ihw{_n}_stratum"]
    _pvalues = [discovery_ihw.loc[_strata == _s, "p_value"].to_numpy() for _s in ihw_strata[_n].index]
    # uniformity above 0.5, five bins per stratum, Bonferroni over the strata
    _tilted = sum(
        stats.chisquare(np.histogram(_p[_p > 0.5], np.linspace(0.5, 1, 6))[0]).pvalue < 0.05 / len(_pvalues)
        for _p in _pvalues
    )
    plot_pvalues_by_stratum(
        _pvalues,
        [f"≤ {_kb:.0f} kb" for _kb in ihw_strata[_n]["to_kb"]],
        title="Above 0.5 the p-values are flat in every stratum"
        if _tilted == 0
        else f"Above 0.5 the p-values are not flat in {_tilted} of {len(_pvalues)} strata",
        subtitle=f"Day 0: pairs per p-value bin, {_n} strata of distance to the TSS; dashed, the count if all were null",
    )
    return


@app.cell
def ihw5_pvalues_plot(
    discovery_ihw,
    ihw_strata,
    np,
    plot_pvalues_by_stratum,
    stats,
):
    _n = 5
    _strata = discovery_ihw[f"ihw{_n}_stratum"]
    _pvalues = [discovery_ihw.loc[_strata == _s, "p_value"].to_numpy() for _s in ihw_strata[_n].index]
    # uniformity above 0.5, five bins per stratum, Bonferroni over the strata
    _tilted = sum(
        stats.chisquare(np.histogram(_p[_p > 0.5], np.linspace(0.5, 1, 6))[0]).pvalue < 0.05 / len(_pvalues)
        for _p in _pvalues
    )
    plot_pvalues_by_stratum(
        _pvalues,
        [f"≤ {_kb:.0f} kb" for _kb in ihw_strata[_n]["to_kb"]],
        title="Above 0.5 the p-values are flat in every stratum"
        if _tilted == 0
        else f"Above 0.5 the p-values are not flat in {_tilted} of {len(_pvalues)} strata",
        subtitle=f"Day 0: pairs per p-value bin, {_n} strata of distance to the TSS; dashed, the count if all were null",
    )
    return


@app.cell
def ihw20_hits_plot(FDR, ihw_strata, plot_hits_by_stratum):
    _n = 20
    plot_hits_by_stratum(
        ihw_strata[_n][["bh_hits", "ihw_hits"]].set_axis(["Benjamini-Hochberg", "IHW"], axis=1),
        [f"{_kb:.0f}" for _kb in ihw_strata[_n]["to_kb"]],
        xlabel="Distance to the TSS, kb: each stratum's upper edge",
        title=f"IHW finds {ihw_strata[_n]['ihw_hits'].sum()} hits, against "
        f"{ihw_strata[_n]['bh_hits'].sum()} with Benjamini-Hochberg",
        subtitle=f"Day 0: hits at FDR {FDR} per stratum of distance to the TSS, {_n} strata",
    )
    return


@app.cell
def ihw5_hits_plot(FDR, ihw_strata, plot_hits_by_stratum):
    _n = 5
    plot_hits_by_stratum(
        ihw_strata[_n][["bh_hits", "ihw_hits"]].set_axis(["Benjamini-Hochberg", "IHW"], axis=1),
        [f"{_kb:.0f}" for _kb in ihw_strata[_n]["to_kb"]],
        xlabel="Distance to the TSS, kb: each stratum's upper edge",
        title=f"IHW finds {ihw_strata[_n]['ihw_hits'].sum()} hits, against "
        f"{ihw_strata[_n]['bh_hits'].sum()} with Benjamini-Hochberg",
        subtitle=f"Day 0: hits at FDR {FDR} per stratum of distance to the TSS, {_n} strata",
    )
    return


@app.cell(hide_code=True)
def md_ihw_lost(mo):
    mo.md(r"""
    #### The hits IHW drops

    IHW drops Benjamini-Hochberg hits where it lowers the weight. Whether those links are
    direct is the question of the long-range section below, and needs its trans pairs.
    """)
    return


@app.cell
def ihw_lost_table(FDR, discovery_ihw, ihw_fits, mo, pd):
    _bh = discovery_ihw["q_value"] < FDR
    _dropped = pd.Series(False, index=discovery_ihw.index)
    _lines = []
    for _n in ihw_fits:
        _lost = _bh & ~(discovery_ihw[f"ihw{_n}_q_value"] < FDR)
        _dropped |= _lost
        _rows = discovery_ihw.loc[_lost]
        _down = _rows["fold_change"] < 1
        _zero = _rows[f"ihw{_n}_weight"] < 1e-6
        _lines.append(
            f"- {_n} strata: IHW drops **{len(_rows)}**, {int(_down.sum())} negative effects and {int((~_down).sum())} "
            f"positive, at a median {_rows['distance_to_tss'].abs().median() / 1000:.0f} kb from the TSS. "
            f"{int(_zero.sum())} have a weight under 10⁻⁶"
            + (f", the smallest p-value among them {_rows.loc[_zero, 'p_value'].min():.1e}." if _zero.any() else ".")
        )
    ihw_lost = discovery_ihw.loc[_dropped].sort_values("distance_to_tss", key=abs)
    mo.vstack(
        [
            mo.md(
                f"Of Benjamini-Hochberg's {int(_bh.sum())} hits:\n\n" + "\n".join(_lines) + "\n\nA weight under 10⁻⁶ "
                "means IHW gave the pair's stratum almost none of the error budget in its fold, so even a very small "
                "p-value does not pass. The table holds the hits dropped with either number of strata."
            ),
            mo.ui.table(ihw_lost, selection=None),
        ]
    )
    return


@app.cell(hide_code=True)
def md_ihw_bins(mo):
    mo.md(r"""
    #### Fewer strata

    Each fold learns its weights from the hits outside it, and 20 strata leave few hits
    per stratum, so a fold can give a stratum no weight at all. Here IHW runs again at 20,
    10, 5 and 3 strata, each with six seeds (the seed sets the folds). Per number of
    strata: the hits, the Benjamini-Hochberg hits lost, the pairs with a weight under
    10⁻⁶, and how many of the lost hits are lost under every seed.
    """)
    return


@app.cell
def ihw_bins(FDR, discovery_ihw, mo, np, pd, pyihw):
    _p = discovery_ihw["p_value"].to_numpy()
    _distance = discovery_ihw["distance_to_tss"].abs().to_numpy()
    _bh = (discovery_ihw["q_value"] < FDR).to_numpy()

    _rows, _lost = [], {}
    for _nbins in (20, 10, 5, 3):
        for _seed in range(1, 7):
            _refit = pyihw.ihw(_p, _distance, FDR, nbins=_nbins, seed=_seed)
            _hit = _refit.adj_pvalues() < FDR
            _lost[_nbins, _seed] = _bh & ~_hit
            _rows.append(
                {
                    "strata": _nbins,
                    "seed": _seed,
                    "ihw_hits": int(_hit.sum()),
                    "gained": int((_hit & ~_bh).sum()),
                    "lost": int(_lost[_nbins, _seed].sum()),
                    "near_zero_pairs": int((_refit.weights() < 1e-6).sum()),
                }
            )
    ihw_bins = pd.DataFrame(_rows)

    _summary = ihw_bins.groupby("strata", sort=False).agg(
        ihw_hits_min=("ihw_hits", "min"),
        ihw_hits_max=("ihw_hits", "max"),
        lost_min=("lost", "min"),
        lost_max=("lost", "max"),
        near_zero_pairs_max=("near_zero_pairs", "max"),
    )
    _every = {_n: np.logical_and.reduce([_lost[_n, _s] for _s in range(1, 7)]) for _n in ihw_bins["strata"].unique()}
    _any = {_n: np.logical_or.reduce([_lost[_n, _s] for _s in range(1, 7)]) for _n in ihw_bins["strata"].unique()}
    _summary["lost_every_seed"] = [int(_every[_n].sum()) for _n in _summary.index]
    _summary["lost_some_seeds"] = [int((_any[_n] & ~_every[_n]).sum()) for _n in _summary.index]
    mo.vstack(
        [
            mo.md(f"Benjamini-Hochberg: **{int(_bh.sum())}** hits. IHW over six seeds per number of strata:"),
            _summary,
            mo.ui.table(ihw_bins, selection=None),
        ]
    )
    return


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
    for _min_cells in sorted(_m for _m in {N_NONZERO_TRT, 10, 15, 20, 30, 50, 100, 200} if _m >= N_NONZERO_TRT):
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
    /// note | Weighting pairs instead of cutting them

    The sweep drops pairs below a threshold, and picking the threshold from its hits
    would undo the FDR control. Independent hypothesis weighting (IHW; Ignatiadis
    et al., 2016) instead gives some pairs more of the error budget than others, with
    no threshold to choose and the FDR still controlled. It runs after discovery, with
    pyIHW, weighted by distance to the TSS. Power, the covariate this sweep uses, is
    not usable there: IHW needs a covariate independent of the p-value under the null,
    and the safe-targeting calibration's p-values run small for highly expressed genes.
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
        response_matrix,
        gene_ids,
        covariate_matrix,
        batch_width=1,
        n_jobs=N_JOBS,
        # n_cells x covariates² floats, built once and shared: each job would otherwise build its own
        x_outer_flat_shared=pysceptre.glm.irls.x_outer_flat(covariate_matrix),
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


@app.cell(hide_code=True)
def md_trans_pairs_qc(mo):
    mo.md(r"""
    ### Which trans pairs can be tested

    Every element in the correction against every gene of its cis pairs that lies on
    another chromosome, chrY left out. A trans pair gets the same pairwise QC as a
    cis pair: a pair too sparse to test in cis would still count toward the trans
    level, and lower it. The cuts start at the cis pairs'; `TRANS_N_NONZERO_TRT`
    and `TRANS_N_NONZERO_CNTRL` change them.
    """)
    return


@app.cell
def trans_pairs_qc(
    N_NONZERO_CNTRL,
    N_NONZERO_TRT,
    discovery_pairs,
    elements,
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
    TRANS_N_NONZERO_TRT, TRANS_N_NONZERO_CNTRL = N_NONZERO_TRT, N_NONZERO_CNTRL  # the cis pairs' cuts

    # the cis pairs in the correction that pass QC give the elements and genes
    _passing = pairs_qc.loc[pairs_qc["pass_qc"], ["grna_target", "response_id"]]
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

    # the cis pairs' QC function
    trans_pairs_qc = pysceptre.pipeline.power.annotate_pairwise_qc(
        _trans_pairs,
        response_matrix,
        gene_ids,
        grna_target_cells,
        response_matrix.shape[1],
        n_nonzero_trt_thresh=TRANS_N_NONZERO_TRT,
        n_nonzero_cntrl_thresh=TRANS_N_NONZERO_CNTRL,
    )
    mo.md(
        f"**{int(trans_pairs_qc['pass_qc'].sum()):,} of {len(trans_pairs_qc):,} trans pairs pass** "
        f"({len(_targets):,} elements by {len(_genes)} genes) at {TRANS_N_NONZERO_TRT} treated and "
        f"{TRANS_N_NONZERO_CNTRL:,} control cells."
    )
    return TRANS_N_NONZERO_CNTRL, TRANS_N_NONZERO_TRT, trans_pairs_qc


@app.cell
def trans_pairs_qc_scatter(
    TRANS_N_NONZERO_CNTRL,
    TRANS_N_NONZERO_TRT,
    mo,
    plot_pairwise_qc,
    trans_pairs_qc,
):
    _fail = int((~trans_pairs_qc["pass_qc"]).sum())
    _on_control = int((trans_pairs_qc["n_nonzero_cntrl"] < TRANS_N_NONZERO_CNTRL).sum())
    mo.vstack(
        [
            plot_pairwise_qc(
                trans_pairs_qc,
                TRANS_N_NONZERO_TRT,
                TRANS_N_NONZERO_CNTRL,
                title=f"{_fail:,} trans pairs fail, "
                + (f"{_on_control:,} of them on the control cut" if _on_control else "all on the treated cut"),
                subtitle="Day 0: trans pairs, cells with a nonzero count of the pair's gene; log scales",
            ),
            mo.md(
                "One dot per trans pair. Treated: cells carrying any guide of the pair's element. Control: all the "
                "other cells. Each axis counts only the cells with a nonzero count of the pair's gene. "
                f"Dashed lines: the cuts, {TRANS_N_NONZERO_TRT} treated and {TRANS_N_NONZERO_CNTRL:,} control cells; "
                "a pair below either line fails."
            ),
        ]
    )
    return


@app.cell
def trans_results(
    N_JOBS,
    SEED,
    SIDE,
    covariate_matrix,
    discovery_pairs,
    gene_fits,
    gene_ids,
    grna_target_cells,
    mo,
    np,
    pairs_qc,
    pd,
    pysceptre,
    response_matrix,
    trans_pairs_qc,
):
    # a submodule pysceptre does not import itself
    from pysceptre.pipeline.trans import trans_fit_chunks as _trans_fit_chunks

    TRANS_CHUNKS_PER_CALL = 4  # discovery's fit chunks per trans call: fewer holds less at once

    # discovery_run's pairs, in its order: the trans test fits each target with the same chunk-mates
    _passing = pairs_qc.loc[pairs_qc["pass_qc"], ["grna_target", "response_id"]]
    _reference = discovery_pairs[["response_id", "grna_target"]].merge(_passing, on=["grna_target", "response_id"])

    # the trans pairs that pass their QC
    _trans_pairs = trans_pairs_qc.loc[trans_pairs_qc["pass_qc"], ["response_id", "grna_target"]].copy()

    # discovery's target fit chunks, rebuilt from its pairs and its settings (run_discovery_analysis's defaults)
    _chunks = _trans_fit_chunks(
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
        f"**{len(trans_results):,}** trans pairs ({_trans_pairs['grna_target'].nunique():,} elements by "
        f"{_trans_pairs['response_id'].nunique()} genes) in "
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
                subtitle="Day 0, Benjamini-Hochberg: pairs significant with a negative effect, by distance to the TSS",
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
def md_corrections_trans(mo):
    mo.md(r"""
    ### Benjamini-Hochberg against IHW, at the trans level

    The same comparison for each correction. IHW calls each pair at its own cutoff, its
    weight times one threshold, so each pair is matched to its gene's trans rate at that
    cutoff rather than at Benjamini-Hochberg's. `above_background_at_cutoffs` does this.
    With one cutoff it reproduces pysceptre's table above, intervals included, and the
    cell checks that. Every element counts, the broad-effect ones included.
    """)
    return


@app.cell
def corrections_trans(
    FDR,
    IHW_STRATA,
    LONG_RANGE_BINS,
    SEED,
    above_background_at_cutoffs,
    discovery_ihw,
    ihw_fits,
    links_vs_trans,
    long_range,
    mo,
    np,
    pd,
    pysceptre,
):
    _keys = ["grna_target", "response_id"]
    _links = long_range.links
    # each cis pair's IHW weights and calls, joined by its element and gene
    _pairs = _links[_keys].merge(
        discovery_ihw[_keys + ["fold_change"] + [f"ihw{_n}_{_c}" for _n in IHW_STRATA for _c in ("weight", "q_value")]],
        on=_keys,
        how="left",
        validate="one_to_one",
    )
    assert _pairs.notna().all().all(), "a cis pair is missing from discovery_ihw"

    # each pair's nominal cutoff: one for Benjamini-Hochberg, its weight times one threshold for IHW
    _cutoffs = {
        "Benjamini-Hochberg": np.full(
            len(_links), pysceptre.analytical_power.bh_nominal_cutoff(discovery_ihw["p_value"].to_numpy(), FDR)
        )
    }
    for _n, _fit in ihw_fits.items():
        _table = _fit.to_dataframe()
        _threshold = _table.loc[_table["adj_pvalue"] < FDR, "weighted_pvalue"].max()  # on p / weight
        _cutoffs[f"IHW, {_n} strata"] = _pairs[f"ihw{_n}_weight"].to_numpy() * _threshold
    _down = _links["log_2_fold_change"].to_numpy() < 0
    _link = {_k: (_links["p_value"].to_numpy() <= _c) & _down for _k, _c in _cutoffs.items()}
    # a correction's links are its negative-effect hits
    assert (_link["Benjamini-Hochberg"] == _links["link"].to_numpy()).all()
    for _n in IHW_STRATA:
        _hits = ((_pairs[f"ihw{_n}_q_value"] < FDR) & (_pairs["fold_change"] < 1)).to_numpy()
        assert (_link[f"IHW, {_n} strata"] == _hits).all(), _n
    # a trans test stopped early counts as no call, exact only while its p-value floor is above every cutoff
    _background = long_range.background
    _early = ((_background["screen"] == "early_stop") & (_background["log_2_fold_change"] < 0)).to_numpy()
    assert (_background.loc[_early, "p_value_floor"] > max(_c.max() for _c in _cutoffs.values())).all()

    # one seed for all: the same bootstrap draws, and for Benjamini-Hochberg pysceptre's own
    corrections_by_distance = {
        _k: above_background_at_cutoffs(
            _links, _background, _cutoffs[_k], _link[_k], distance_bins=LONG_RANGE_BINS, seed=SEED
        )
        for _k in _cutoffs
    }
    pd.testing.assert_frame_equal(
        corrections_by_distance["Benjamini-Hochberg"], long_range.by_distance, check_dtype=False, rtol=1e-9
    )

    _counts = {_k: links_vs_trans(_t) for _k, _t in corrections_by_distance.items()}
    _totals = pd.DataFrame(
        {_k: _c[["significant cis", "significant trans, expected", "cis above the trans level"]].sum() for _k, _c in _counts.items()}
    ).T.round(1)
    _bins = list(_counts["Benjamini-Hochberg"].index)
    _by_bin = pd.concat(_counts, names=["correction"]).swaplevel()
    _by_bin = _by_bin.loc[[(_b, _k) for _b in _bins for _k in _counts]]

    mo.vstack(
        [
            mo.md("Totals over the distance bins: links, the trans level's expected links, and the links above it."),
            _totals,
            _by_bin,
        ]
    )
    return (corrections_by_distance,)


@app.cell
def ihw20_long_range_plot(
    LONG_RANGE_BINS,
    corrections_by_distance,
    links_vs_trans,
    np,
    plot_links_vs_trans,
):
    _n = 20
    _counts = links_vs_trans(corrections_by_distance[f"IHW, {_n} strata"])
    # the bins within 100 kb against those beyond it
    _upper = np.asarray(LONG_RANGE_BINS[1:], dtype=float)
    _near, _far = _counts[_upper <= 100_000], _counts[_upper > 100_000]


    def _share(d):
        _k = d["significant cis"].sum()
        return f"{d['cis above the trans level'].sum() / _k:.0%}" if _k else "none"


    plot_links_vs_trans(
        _counts,
        title=f"{_share(_near)} of links within 100 kb clear the trans level, {_share(_far)} beyond",
        subtitle=f"Day 0, IHW with {_n} strata: pairs significant with a negative effect, by distance to the TSS",
    )
    return


@app.cell
def ihw5_long_range_plot(
    LONG_RANGE_BINS,
    corrections_by_distance,
    links_vs_trans,
    np,
    plot_links_vs_trans,
):
    _n = 5
    _counts = links_vs_trans(corrections_by_distance[f"IHW, {_n} strata"])
    # the bins within 100 kb against those beyond it
    _upper = np.asarray(LONG_RANGE_BINS[1:], dtype=float)
    _near, _far = _counts[_upper <= 100_000], _counts[_upper > 100_000]


    def _share(d):
        _k = d["significant cis"].sum()
        return f"{d['cis above the trans level'].sum() / _k:.0%}" if _k else "none"


    plot_links_vs_trans(
        _counts,
        title=f"{_share(_near)} of links within 100 kb clear the trans level, {_share(_far)} beyond",
        subtitle=f"Day 0, IHW with {_n} strata: pairs significant with a negative effect, by distance to the TSS",
    )
    return


@app.cell(hide_code=True)
def md_direct_indirect(mo):
    mo.md(r"""
    ### Direct and indirect effects, as in DC-TAP

    The DC-TAP paper gives each link a probability of being a direct effect, from its distance to the
    TSS alone (the lab's `CRISPR_indirect_effects` workflow):

    1. **Indirect rate**, one number: the share of trans tests (elements against genes on other
       chromosomes) significant with a negative effect, each called at the cis cutoff, the largest
       p-value Benjamini-Hochberg calls significant.
    2. **Cis rate**: the same share among cis pairs, in 50 kb bins of distance to the TSS.
    3. **Direct rate** of a bin: its cis rate minus the indirect rate, floored at 0. A power law,
       log(direct rate) against log(mean distance), is fitted to the bins above 0.
    4. A link at distance *d* is direct with probability direct(*d*) / (direct(*d*) + indirect).

    As in DC-TAP, the rates come from distal elements, here those neither a promoter nor a labelled
    enhancer, and the probability is then read for every link. It is DC-TAP's model for negative
    effects; positive-effect hits are left out. Three differences: DC-TAP averaged the direct rates of
    six screens before fitting, and this is one screen; DC-TAP tested 100 random trans genes per
    element, and here every trans test pysceptre ran counts; DC-TAP kept the genes of its TSS file,
    which is not applied here.

    The probability depends on distance alone, and so do IHW's weights, so the links IHW drops sit at
    low probability by construction. What the cell shows is how low.
    """)
    return


@app.cell
def direct_indirect(
    FDR,
    IHW_STRATA,
    discovery_ihw,
    discovery_table,
    long_range,
    mo,
    np,
    pd,
    pysceptre,
):
    # DC-TAP's direct/indirect model (EngreitzLab/CRISPR_indirect_effects), on Benjamini-Hochberg's calls
    _keys = ["grna_target", "response_id"]
    _links = long_range.links

    # the rates come from distal elements: neither a promoter nor a labelled enhancer
    _types = discovery_table.loc[
        discovery_table["grna_target"].isin(_links["grna_target"]), ["grna_target", "genomic_element"]
    ].drop_duplicates()
    assert not _types["grna_target"].duplicated().any(), "an element with more than one genomic_element"
    assert set(_types["grna_target"]) == set(_links["grna_target"])
    _distal = set(_types.loc[~_types["genomic_element"].isin(["promoter", "enhancer"]), "grna_target"])


    def _called(d, c):
        """A call at cutoff c: p <= c with a negative effect. A NaN p (screened out) is tested, not called."""
        return ((d["p_value"] <= c) & (d["log_2_fold_change"] < 0)).to_numpy()


    # DC-TAP's cutoff: the largest p-value Benjamini-Hochberg calls significant, in either direction
    _p = discovery_ihw["p_value"].to_numpy()
    _cutoff = _p[(discovery_ihw["q_value"] < FDR).to_numpy()].max()
    assert (_called(_links, _cutoff) == _links["link"].to_numpy()).all()
    _bh_cutoff = pysceptre.analytical_power.bh_nominal_cutoff(_p, FDR)
    _trans_calls = {_c: int(_called(long_range.background, _c).sum()) for _c in (_cutoff, _bh_cutoff)}
    _trans = long_range.background[long_range.background["grna_target"].isin(_distal)]
    _indirect = float(_called(_trans, _cutoff).mean())

    # cis rates in 50 kb bins, and the power law through the bins whose direct rate is above 0
    _cis = _links[_links["grna_target"].isin(_distal)]
    _d = _cis["distance"].abs()
    _edges = np.arange(0, 1_000_001, 50_000)
    _bin = pd.cut(
        _d, _edges, include_lowest=True, labels=[f"{a / 1e3:g} to {b / 1e3:g} kb" for a, b in zip(_edges[:-1], _edges[1:])]
    )
    assert _bin.notna().all(), "a distal pair beyond 1 Mb"
    _rows = _cis.groupby(_bin, observed=False)
    dctap_bins = pd.DataFrame(
        {"pairs": _rows.size(), "links": _rows["link"].sum(), "mean_distance": _d.groupby(_bin, observed=False).mean()}
    )
    dctap_bins["direct_rate"] = (dctap_bins["links"] / dctap_bins["pairs"] - _indirect).clip(lower=0)
    _fitted = dctap_bins["direct_rate"] > 0
    _b, _log_a = np.polyfit(np.log(dctap_bins.loc[_fitted, "mean_distance"]), np.log(dctap_bins.loc[_fitted, "direct_rate"]), 1)
    dctap_fit = {"a": float(np.exp(_log_a)), "b": float(_b), "indirect_rate": _indirect, "cutoff": float(_cutoff)}


    def _direct_rate(distance):
        return np.minimum(dctap_fit["a"] * np.abs(np.asarray(distance, dtype=float)) ** dctap_fit["b"], 1)


    def _p_direct(distance):
        return _direct_rate(distance) / (_direct_rate(distance) + _indirect)


    dctap_bins["indirect_expected"] = dctap_bins["pairs"] * _indirect
    dctap_bins["direct_fitted"] = pd.Series(_direct_rate(_d), index=_cis.index).groupby(_bin, observed=False).sum()
    _grid = np.arange(1_000, 1_000_001, 1_000)
    dctap_curve = pd.DataFrame({"distance": _grid, "p_direct": _p_direct(_grid)})

    # every Benjamini-Hochberg link, any element, with its probability and what each IHW does with it
    _hits = _links.loc[_links["link"], _keys + ["distance"]].merge(
        discovery_ihw[_keys + ["gene_symbol", "genomic_element"] + [f"ihw{_n}_q_value" for _n in IHW_STRATA]],
        on=_keys,
        how="left",
        validate="one_to_one",
    )
    _drop20, _drop5 = (_hits["ihw20_q_value"] >= FDR).to_numpy(), (_hits["ihw5_q_value"] >= FDR).to_numpy()
    DCTAP_IHW = ["kept by both", "dropped by IHW 20 only", "dropped by IHW 5 only", "dropped by both"]
    dctap_links = _hits.assign(
        distance_kb=_hits["distance"].abs() / 1e3,
        p_direct=_p_direct(_hits["distance"]),
        ihw=np.select(
            [_drop20 & _drop5, _drop20, _drop5],
            ["dropped by both", "dropped by IHW 20 only", "dropped by IHW 5 only"],
            default="kept by both",
        ),
    )

    # the fit against the counts it was fitted to: the links above the indirect level, all and beyond 500 kb
    _above = dctap_bins["links"] - dctap_bins["indirect_expected"]
    _far = (dctap_bins["mean_distance"] > 500_000).to_numpy()
    _check = pd.DataFrame(
        {
            "links above the indirect level": [_above.sum(), _above[_far].sum()],
            "direct links the fit expects": [dctap_bins["direct_fitted"].sum(), dctap_bins.loc[_far, "direct_fitted"].sum()],
        },
        index=["distal pairs, all", "distal pairs beyond 500 kb"],
    ).round(1)
    _summary = (
        dctap_links.groupby("ihw")
        .agg(
            links=("p_direct", "size"),
            median_kb=("distance_kb", "median"),
            median_p_direct=("p_direct", "median"),
            p_direct_below_half=("p_direct", lambda _x: int((_x < 0.5).sum())),
        )
        .reindex(DCTAP_IHW)
        .round(2)
    )
    _above90 = dctap_curve.loc[dctap_curve["p_direct"] >= 0.9, "distance"]
    _above50 = dctap_curve.loc[dctap_curve["p_direct"] >= 0.5, "distance"]
    mo.vstack(
        [
            mo.md(
                f"Indirect rate **{_indirect:.3%}**: {int(_called(_trans, _cutoff).sum())} of {len(_trans):,} trans tests of "
                f"{len(_distal):,} distal elements, called at p ≤ {_cutoff:.2e}"
                + (
                    ""
                    if _trans_calls[_cutoff] == _trans_calls[_bh_cutoff]
                    else f" (pysceptre's cutoff, {_bh_cutoff:.2e}, calls {_trans_calls[_bh_cutoff]} trans tests to DC-TAP's {_trans_calls[_cutoff]})"
                )
                + f". {int(_fitted.sum())} of {len(dctap_bins)} bins have a direct rate above 0 and enter the fit: "
                f"direct rate = {dctap_fit['a']:.3g} × distance^{dctap_fit['b']:.2f}. "
                f"A link is at least 90% likely direct up to "
                f"{(f'{_above90.max() / 1e3:g} kb' if len(_above90) else 'no distance')}, and at least 50% up to "
                f"{(f'{_above50.max() / 1e3:g} kb' if len(_above50) else 'no distance')}."
            ),
            dctap_bins.drop(columns="mean_distance").round({"direct_rate": 5, "indirect_expected": 2, "direct_fitted": 2}),
            mo.md("The fit against the counts it was fitted to:"),
            _check,
            mo.md("Benjamini-Hochberg's links by what IHW does with them:"),
            _summary,
            dctap_links.loc[dctap_links["ihw"] != DCTAP_IHW[0], ["gene_symbol", "genomic_element", "distance_kb", "p_direct", "ihw"]]
            .sort_values("distance_kb")
            .round({"distance_kb": 1, "p_direct": 2})
            .reset_index(drop=True),
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


    return headline, house_style


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


@app.cell(hide_code=True)
def ihw_plot_fn(np, plt):
    def plot_ihw_weights(weights, labels, xlabel="", title="", subtitle=""):
        """IHW's weight per stratum: each fold a thin grey line, their mean in blue, 1 (Benjamini-Hochberg's) in black.

        `weights`: strata x folds, as `IHWResult.weights(levels_only=True)` returns it.
        `labels`: one tick label per stratum.
        """
        strata = np.arange(len(weights))
        height = 4.4
        fig, ax = plt.subplots(figsize=(7.5, height))
        fig.subplots_adjust(left=0.08, right=0.97, bottom=0.17, top=1 - 1.05 / height)
        ax.axhline(1, color="0.1", linewidth=1, label="Benjamini-Hochberg, 1")
        for fold, column in enumerate(weights.T):
            ax.plot(strata, column, color="0.7", linewidth=0.8, marker="o", markersize=2.5,
                    label="One fold" if fold == 0 else None)
        ax.plot(strata, weights.mean(axis=1), color="#0072B2", linewidth=2, marker="o", markersize=4,
                label="Mean of the folds")
        ax.set_xticks(strata, labels)
        ax.set_xlim(-0.6, len(strata) - 0.4)
        ax.set_ylim(0, None)
        ax.set_xlabel(xlabel, fontsize=9)
        ax.set_axisbelow(True)
        ax.grid(axis="y", linewidth=0.5, color="#B7C6CF")
        ax.grid(axis="x", visible=False)
        ax.spines[["top", "right", "left", "bottom"]].set_visible(False)
        ax.tick_params(axis="both", length=0, labelsize=8.5)
        ax.legend(loc="lower left", bbox_to_anchor=(0, 1.02), ncol=3, frameon=False, fontsize=8.5,
                  handlelength=1.6, borderaxespad=0, borderpad=0)
        fig.text(0.03, 1 - 0.3 / height, title, fontsize=13, fontweight="bold", ha="left")
        fig.text(0.03, 1 - 0.55 / height, subtitle, fontsize=9.5, color="0.3", ha="left")
        return fig


    def plot_pvalues_by_stratum(pvalues, labels, bins=10, ncols=5, title="", subtitle=""):
        """P-value counts per stratum in small multiples, each with the count a bin holds if every pair is null.

        `pvalues`: one array per stratum.
        `labels`: one panel title per stratum.
        """
        edges = np.linspace(0, 1, bins + 1)
        counts = [np.histogram(p, edges)[0] for p in pvalues]
        nrows = -(-len(pvalues) // ncols)
        height = 1.3 + 1.25 * nrows
        fig, axes = plt.subplots(nrows, ncols, figsize=(7.5, height), sharex=True, sharey=True, squeeze=False)
        fig.subplots_adjust(left=0.08, right=0.97, bottom=0.45 / height, top=1 - 1.05 / height, hspace=0.55, wspace=0.15)
        top = max(c.max() for c in counts) * 1.08
        for ax, count, p, label in zip(axes.flat, counts, pvalues, labels):
            ax.bar(edges[:-1], count, width=1 / bins, align="edge", color="#0072B2", linewidth=0)
            ax.axhline(len(p) / bins, color="0.1", linewidth=1, linestyle=(0, (3, 2)))
            ax.set_title(label, fontsize=8.5, loc="left", pad=3)
            ax.set_xlim(0, 1)
            ax.set_ylim(0, top)
            ax.set_xticks([0, 0.5, 1], ["0", "0.5", "1"])
            ax.set_axisbelow(True)
            ax.grid(axis="y", linewidth=0.5, color="#B7C6CF")
            ax.grid(axis="x", visible=False)
            ax.spines[["top", "right", "left", "bottom"]].set_visible(False)
            ax.tick_params(axis="both", length=0, labelsize=8.5)
        for ax in axes.flat[len(pvalues):]:
            ax.set_visible(False)
        fig.text(0.03, 1 - 0.3 / height, title, fontsize=13, fontweight="bold", ha="left")
        fig.text(0.03, 1 - 0.55 / height, subtitle, fontsize=9.5, color="0.3", ha="left")
        return fig


    def plot_hits_by_stratum(hits, labels, xlabel="", title="", subtitle=""):
        """Hits per stratum, one bar per correction side by side, each numbered.

        `hits`: one row per stratum, one column per correction (at most two), counts.
        `labels`: one tick label per stratum.
        """
        # Okabe-Ito, both at least 3:1 on white
        colours = ["#D55E00", "#0072B2"]
        if hits.shape[1] > len(colours):
            raise ValueError(f"At most {len(colours)} corrections, not {hits.shape[1]}.")
        strata = np.arange(len(hits))
        width = 0.4
        height = 4.4
        fig, ax = plt.subplots(figsize=(7.5, height))
        fig.subplots_adjust(left=0.08, right=0.97, bottom=0.17, top=1 - 1.05 / height)
        for i, ((name, count), colour) in enumerate(zip(hits.items(), colours)):
            x = strata + (i - 0.5) * width
            ax.bar(x, count, width=width, color=colour, linewidth=0, label=name)
            for xi, n in zip(x, count):
                if n:
                    ax.annotate(f"{n:.0f}", (xi, n), xytext=(0, 2), textcoords="offset points", ha="center",
                                va="bottom", fontsize=7, color="0.1")
        ax.set_xticks(strata, labels)
        ax.set_xlim(-0.6, len(strata) - 0.4)
        ax.set_ylim(0, hits.to_numpy().max() * 1.1)
        ax.set_xlabel(xlabel, fontsize=9)
        ax.set_axisbelow(True)
        ax.grid(axis="y", linewidth=0.5, color="#B7C6CF")
        ax.grid(axis="x", visible=False)
        ax.spines[["top", "right", "left", "bottom"]].set_visible(False)
        ax.tick_params(axis="both", length=0, labelsize=8.5)
        ax.legend(loc="lower left", bbox_to_anchor=(0, 1.02), ncol=hits.shape[1], frameon=False, fontsize=8.5,
                  handlelength=1.2, borderaxespad=0, borderpad=0)
        fig.text(0.03, 1 - 0.3 / height, title, fontsize=13, fontweight="bold", ha="left")
        fig.text(0.03, 1 - 0.55 / height, subtitle, fontsize=9.5, color="0.3", ha="left")
        return fig


    return plot_hits_by_stratum, plot_ihw_weights, plot_pvalues_by_stratum


@app.cell(hide_code=True)
def trans_level_fn(np, pd, sparse):
    def above_background_at_cutoffs(links, background, cutoff, link, *, distance_bins, n_bootstrap=2000, seed=None):
        """pysceptre's `above_background_by_distance`, with each cis pair called at its own cutoff.

        IHW calls each pair at its weight times one threshold. A pair's matched background is then its gene's
        trans call rate at that pair's cutoff. With one cutoff for every pair this is pysceptre's function, draws
        included.

        `links`, `background`: `run_specificity_check`'s tables of the same names.
        `cutoff`, `link`: per row of `links`, the nominal p-value it is called at, and whether it is a link.

        Returns pysceptre's table: per bin `distance_bin`, `pairs`, `links`, then `cis_rate`, `background` and
        `above`, each with `_lo` and `_hi` (the 2.5 and 97.5 percentiles of a bootstrap over elements).
        """
        import warnings

        edges = np.asarray(distance_bins, dtype=float)
        n_bins = edges.size - 1
        elements = pd.Index(np.sort(links["grna_target"].unique()))
        genes = pd.Index(np.sort(links["response_id"].unique()))
        shape = (len(elements), len(genes))
        rng = np.random.default_rng(seed)
        weights = np.vstack(
            [np.ones(len(elements)), rng.multinomial(len(elements), np.full(len(elements), 1.0 / len(elements)), n_bootstrap)]
        ).astype(float)

        def per_gene(rows, cols):
            """Each draw's weighted count per gene of the (element, gene) entries given."""
            return (sparse.csc_matrix((np.ones(rows.size), (rows, cols)), shape=shape).T @ weights.T).T

        te = elements.get_indexer(background["grna_target"])
        tg = genes.get_indexer(background["response_id"])
        inside = (te >= 0) & (tg >= 0)
        te, tg = te[inside], tg[inside]
        p = background["p_value"].to_numpy(dtype=float)[inside]
        down = background["log_2_fold_change"].to_numpy(dtype=float)[inside] < 0
        tested = per_gene(te, tg)

        ce = elements.get_indexer(links["grna_target"])
        cg = genes.get_indexer(links["response_id"])
        cb = pd.cut(links["distance"], edges, include_lowest=True, labels=False).to_numpy()
        if np.isnan(cb).any():
            raise ValueError(f"{int(np.isnan(cb).sum())} pair(s) have a distance outside distance_bins")
        cb = cb.astype(int)
        pairs_eb, links_eb = np.zeros((len(elements), n_bins)), np.zeros((len(elements), n_bins))
        np.add.at(pairs_eb, (ce, cb), 1.0)
        np.add.at(links_eb, (ce, cb), np.asarray(link, dtype=float))

        cutoff = np.asarray(cutoff, dtype=float)
        expected = np.zeros((weights.shape[0], n_bins))
        for c in np.unique(cutoff):
            # a trans test is a call at c as a cis pair is: p <= c with a negative effect
            called = down & (p <= c)
            with np.errstate(divide="ignore", invalid="ignore"):
                rate = np.nan_to_num(per_gene(te[called], tg[called]) / tested)
            rows = np.flatnonzero(cutoff == c)
            for b in range(n_bins):
                r = rows[cb[rows] == b]
                expected[:, b] += (per_gene(ce[r], cg[r]) * rate).sum(axis=1)

        n = weights @ pairs_eb
        out = pd.DataFrame(
            {
                "distance_bin": [
                    f"> {a / 1000:g} kb" if np.isinf(b) else f"{a / 1000:g} to {b / 1000:g} kb"
                    for a, b in zip(edges[:-1], edges[1:])
                ],
                "pairs": pairs_eb.sum(0).astype(int),
                "links": links_eb.sum(0).astype(int),
            }
        )
        with np.errstate(divide="ignore", invalid="ignore"), warnings.catch_warnings():
            # a bin no draw can estimate is NaN, as in pysceptre
            warnings.simplefilter("ignore", RuntimeWarning)
            cis, back = (weights @ links_eb) / n, expected / n
            for name, d in (("cis_rate", cis), ("background", back), ("above", (cis - back) / cis)):
                out[name] = d[0]
                out[f"{name}_lo"], out[f"{name}_hi"] = np.nanpercentile(d[1:], [2.5, 97.5], axis=0)
        return out


    return (above_background_at_cutoffs,)


if __name__ == "__main__":
    app.run()
