import marimo

__generated_with = "0.19.2"
app = marimo.App(width="full")


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell
def _():
    from datetime import datetime
    from pathlib import Path
    from pyvis.network import Network
    from statsmodels.stats.multitest import multipletests

    import altair as alt

    import pandas as pd
    import networkx as nx
    import matplotlib

    def fmt_int(n):
        return f"{int(n):,}".replace(",", "'")

    def fmt_float(x):
        return "—" if pd.isna(x) else f"{x:.2f}"
    return Network, Path, alt, datetime, fmt_float, fmt_int, multipletests, pd


@app.cell
def _(datetime, fmt_float, fmt_int, mo, results_flt, summary):
    _left = mo.md(
        f"""
    **Date:** {datetime.strptime("2025/12/09", "%Y/%m/%d").strftime("%B %d, %Y")}  
    **Project:** hiPSC -> EC differentiaton    
    **Timepoints:** 3
    """
    )

    _right = mo.md(
        f"""
    **Total pairs tested:** {fmt_int(results_flt.shape[0])}  
    **Total significant pairs:** {fmt_int(results_flt["significant"].value_counts().get(True, 0))}  
    **Tested genes:** {fmt_int(results_flt["response_id"].nunique())}  
    """
    ).right()

    _header = mo.md(
        f"""
    ## Results Overview
    {mo.hstack([_left, _right], widths="equal")}
    """
    )

    rows = []
    # 
    for _, row_summary in summary.sort_values("day").iterrows():
        day = row_summary["day"]

        day_cards = [
            mo.stat(
                label=f"Day {day} · Total significant pairs",
                value=fmt_int(row_summary["n_significant"]),
                bordered=True,
            ).style(min_width="240px"),

            mo.stat(
                label=f"Day {day} · significant",
                value=fmt_int(row_summary["n_neg"]),
                bordered=True,
                direction="decrease",
                caption="Negative effect",
            ).style(min_width="240px"),

            mo.stat(
                label=f"Day {day} · mean effect",
                value=fmt_float(row_summary["mean_effect_neg"]),
                bordered=True,
                direction="decrease",
                caption="Negative effect",
            ).style(min_width="240px"),

            mo.stat(
                label=f"Day {day} · unique genes",
                value=fmt_int(row_summary["unique_gene_neg"]),
                bordered=True,
                direction="decrease",
                caption="Decreased transcription",
            ).style(min_width="240px"),

            mo.stat(
                label=f"Day {day} · unique elements",
                value=fmt_int(row_summary["unique_element_neg"]),
                bordered=True,
                direction="decrease",
                caption="Negative effect",
            ).style(min_width="240px"),

            mo.stat(
                label=f"Day {day} · significant",
                value=fmt_int(row_summary["n_pos"]),
                bordered=True,
                direction="increase",
                caption="Positive effect",
            ).style(min_width="240px"),

            mo.stat(
                label=f"Day {day} · mean effect",
                value=fmt_float(row_summary["mean_effect_pos"]),
                bordered=True,
                direction="increase",
                caption="Positive effect",
            ).style(min_width="240px"),

            mo.stat(
                label=f"Day {day} · unique genes",
                value=fmt_int(row_summary["unique_gene_pos"]),
                bordered=True,
                direction="increase",
                caption="Increased transcription",
            ).style(min_width="240px"),

            mo.stat(
                label=f"Day {day} · unique elements",
                value=fmt_int(row_summary["unique_element_pos"]),
                bordered=True,
                direction="increase",
                caption="Positive effect",
            ).style(min_width="240px"),
        ]

        rows.append(
            mo.hstack(
                day_cards,
                align="center",
                justify="space-around",
            )
        )

    mo.vstack([_header, *rows])
    return


@app.cell
def _(alt, mo, results_flt):
    # For each (reponse_id, grna_target) pair, count how many days it was significant
    sig_counts = (
        results_flt
        .groupby(["response_id", "grna_target"])
        .agg(
            n_significant=("significant", "sum"),
            total_days=("day", "nunique")
        )
        .reset_index()
    )
    sig_days = (
        results_flt[results_flt["significant"]]
        .groupby(["response_id", "grna_target"])["day"]
        .agg(lambda x: ",".join(sorted(x.astype(str).unique())))
        .reset_index(name="sig_days")
    )

    # 4) Combine both
    sig_counts = sig_counts.merge(sig_days, on=["response_id", "grna_target"], how="left")
    sig_counts.sig_days.value_counts().reset_index()

    bars = (
        alt.Chart(sig_counts.sig_days.value_counts().reset_index())
        .mark_bar()
        .encode(
            x=alt.X("count:Q", title="Number of significant pairs"),
            y=alt.Y("sig_days:N", sort="-x", title="Intersection of days")
        )
        .properties(
            width=500,
            title={
                "text": "Number of significant (response_id, grna_target) pairs by intersection of days",
                "subtitle": f"Number of unique significant (response_id, grna_target) pairs (n={sig_counts.sig_days.value_counts().sum()})",
                "anchor": "start"
            }
        )
        .configure_axis(
            labelFontSize=14,   # axis tick labels
            titleFontSize=16    # axis titles
        ).add_params(alt.selection_point(name="set"))
        .configure_title(fontSize=16, anchor="start")

    )
    chart_bar = mo.ui.altair_chart(bars)
    chart_bar
    return chart_bar, sig_counts, sig_days


@app.cell
def _(sig_days):
    sig_days
    return


@app.cell
def _(chart_bar, df, mo, sig_counts):
    mo.stop(chart_bar.value is None)
    #sig_days[chart_bar.value["sig_days"]]
    selected_intersection = chart_bar.value["sig_days"].tolist()
    intersection_df = sig_counts[sig_counts["sig_days"].isin(selected_intersection)]
    selected_intersection_table = mo.ui.table(intersection_df, selection="single")

    dropdown_gene = mo.ui.dropdown.from_series(intersection_df["response_id"])
    dropdown_element = mo.ui.dropdown.from_series(df["grna_target"])

    mo.vstack([dropdown_gene,
                  dropdown_element,
                  selected_intersection_table]
                )
    return (selected_intersection_table,)


@app.cell
def _(mo, selected_intersection_table):
    mo.stop(selected_intersection_table is None)
    selected_intersection_table.value
    return


@app.cell
def _(Network, mo, selected_intersection_table, sig_counts_flt):
    # Stop if nothing selected
    v = selected_intersection_table.value
    mo.stop(v is None or v.empty)

    # Create a Network object
    net = Network("800px", "800px", notebook=True, cdn_resources="remote")

    resp_nodes = set()
    grna_nodes = set()

    edge_color_map = {1: "lightgray", 2: "orange", 3: "red"}

    selected_gene = v["response_id"].iloc[0]
    selected_element = v["grna_target"].iloc[0]

    # Filter df
    selected_mask = (sig_counts_flt.response_id == selected_gene) | (sig_counts_flt.grna_target == selected_element)

    for _, row in sig_counts_flt[selected_mask].iterrows():
        resp = row["response_id"]
        grna = row["grna_target"]
        n_days = int(row["n_significant"])
        sig_days_list = row["sig_days"]

        resp_id = f"resp::{resp}"
        grna_id = f"grna::{grna}"

        # (Optional) you no longer need this, since mask already enforces it:
        # if resp != selected_gene and grna != selected_element:
        #     continue

        if resp_id not in resp_nodes:
            net.add_node(
                resp_id,
                label=str(resp),
                color="blue",
                shape="triangle" if resp == selected_gene else "dot",
                title=f"response_id: {resp}",
            )
            resp_nodes.add(resp_id)

        if grna_id not in grna_nodes:
            net.add_node(
                grna_id,
                label=str(grna),
                color="yellow",
                shape="triangle" if grna == selected_element else "dot",
                title=f"gRNA target: {grna}",
            )
            grna_nodes.add(grna_id)

        net.add_edge(
            resp_id,
            grna_id,
            value=n_days,
            color=edge_color_map.get(n_days, "black"),
            title=f"{n_days} significant day(s): {sig_days_list}",
        )

    net.show_buttons(filter_=["physics", "layout", "interaction"])
    html = net.generate_html()

    legend_md = """
    **Nodes:** 
    <span style="color:blue">●</span> Gene &nbsp;
    <span style="color:blue">▲</span> Selected Gene &nbsp;
    <span style="color:gold">●</span> Element &nbsp;
    <span style="color:gold">▲</span> Selected element

    **Edges:** 
    <span style="color:lightgray">━━</span> 1 timepoint &nbsp;
    <span style="color:orange">━━</span> 2 timepoints &nbsp;
    <span style="color:red">━━</span> 3 timepoints
    """

    mo.vstack([mo.md(legend_md), mo.iframe(html)])

    return


@app.cell
def _(results_flt):
    results_flt[results_flt["significant"]].groupby(["grna_target","day"]).size()
    return


@app.cell
def _(sig_counts):
    sig_counts_flt = sig_counts[sig_counts["n_significant"] > 0]
    sig_counts_flt.head()
    return (sig_counts_flt,)


@app.cell
def _(results_flt):
    results_flt[results_flt["significant"]].groupby(["response_id","day"]).size()
    return


@app.cell
def _(results_flt):
    results_flt[results_flt["significant"]].groupby(["grna_target","day"]).size()
    return


@app.cell
def _():
    import numpy as np
    import matplotlib.pyplot as plt

    # ------------------------------------------------------------
    # 1. Create a toy contact matrix with TAD-like structure
    # ------------------------------------------------------------
    np.random.seed(0)

    n = 200  # matrix size
    mat = np.random.poisson(lam=2, size=(n, n)).astype(float)

    # Add block structures (fake TADs)
    for start in [20, 70, 130]:
        end = start + 40
        mat[start:end, start:end] += np.random.poisson(lam=15, size=(40, 40))

    # Symmetrize
    mat = (mat + mat.T) / 2

    # ------------------------------------------------------------
    # 2. Convert to triangular Hi-C map (upper triangle)
    # ------------------------------------------------------------
    tri = np.triu(mat)

    # ------------------------------------------------------------
    # 3. Plot triangular Hi-C contact map
    # ------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(6, 6))

    im = ax.imshow(
        np.log1p(tri),       # log transform improves visibility
        cmap="YlGnBu",       # similar to many Hi-C papers
        origin="lower",
        interpolation="nearest"
    )

    # Remove ticks for a cleaner look
    ax.set_xticks([])
    ax.set_yticks([])

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("log1p(contacts)", fontsize=12)

    ax.set_title("Triangular Hi-C–Style Contact Map", fontsize=14)

    plt.tight_layout()
    plt.show()
    return (np,)


@app.cell
def _(all_days):
    # Keep only the pairs that pass quality control
    results_flt = all_days.loc[list(all_days["pass_qc"]),]
    results_flt.shape, results_flt.head()
    return (results_flt,)


@app.cell
def _():
    return


@app.cell
def _(results_flt):
    summary = (
        results_flt
        .assign(
            # keep effect size only if significant
            effect_size_sig=lambda d: d["effect_size"].where(d["significant"]),

            # split by sign
            effect_pos=lambda d: d["effect_size"].where(d["significant"] & (d["effect_size"] > 0)),
            effect_neg=lambda d: d["effect_size"].where(d["significant"] & (d["effect_size"] < 0)),

            # genes by sign
            gene_pos=lambda d: d["response_id"].where(d["significant"] & (d["effect_size"] > 0)),
            gene_neg=lambda d: d["response_id"].where(d["significant"] & (d["effect_size"] < 0)),

            # element by sign
            element_pos=lambda d: d["grna_target"].where(d["significant"] & (d["effect_size"] > 0)),
            element_neg=lambda d: d["grna_target"].where(d["significant"] & (d["effect_size"] < 0)),
        )
        .groupby("day")
        .agg(
            n_significant=("significant", "sum"),

            n_pos=("effect_pos", "count"),
            n_neg=("effect_neg", "count"),

            mean_effect_pos=("effect_pos", "mean"),
            mean_effect_neg=("effect_neg", "mean"),

            unique_gene_pos=("gene_pos", "nunique"),
            unique_gene_neg=("gene_neg", "nunique"),

            unique_element_pos=("element_pos", "nunique"),
            unique_element_neg=("element_neg", "nunique"),
        )
        .reset_index()
    )
    summary
    return (summary,)


@app.cell
def _(Path, multipletests, np, pd):
    # Load the SCEPTRE results
    base = Path("results/04-2025-10-27_sceptre_analysis/2025-11-13_sceptre_outputs")
    days = ["0", "2", "4"]

    dfs = []
    for d in days:
        df = pd.read_csv(
            base / f"day{d}" / "sceptre_outputs" / "sceptre_discovery_results.csv"
        )
        df["day"] = d

        # FDR on the p_value column for this day
        mask = df["p_value"].notna()
        pvals = df.loc[mask, "p_value"].values

        if pvals.size > 0:
            _, pval_fdr, _, _ = multipletests(pvals, method="fdr_bh")
            df["pvalue_fdr"] = np.nan
            df.loc[mask, "pvalue_fdr"] = pval_fdr
        else:
            df["pvalue_fdr"] = np.nan

        dfs.append(df)

    all_days = pd.concat(dfs, ignore_index=True)
    # Make sure true/false strings are interpreted as boolean
    all_days["significant"] = (
        all_days["significant"]
        .replace({'True': True, 'False': False})   # handles stray strings
        .astype(bool)
    )
    # Make sure true/false strings are interpreted as boolean
    all_days["pass_qc"] = (
        all_days["pass_qc"]
        .replace({'True': True, 'False': False})   # handles stray strings
        .astype(bool)
    )
    # Add effect sizes
    all_days["effect_size"] = (all_days["fold_change"] - 1) * 100
    return all_days, df


@app.cell(disabled=True, hide_code=True)
def _(alt, results_flt):
    mask_not_na = results_flt["p_value"].notna()
    df_plot = results_flt.loc[mask_not_na, ["p_value","pvalue_fdr"]].copy()

    n_sig = (df_plot["pvalue_fdr"] <= 0.1).sum()

    # Plot p-values
    pval_chart = alt.Chart(df_plot).mark_bar().encode(
        x=alt.X("p_value:Q", bin=alt.Bin(maxbins=50), title="Raw p-value"),
        y=alt.Y("count():Q", title="Count")
    ).properties(
        width=500,
        height=200,
        title={
            "text": "Distribution of raw p-values",
            "subtitle": f"Significant hits at 1% FDR highlighted (n={n_sig})",
            "anchor": "start",
        }
    )

    # Plot FDR
    fdr_chart = alt.Chart(df_plot).mark_bar().encode(
        x=alt.X("pvalue_fdr:Q", bin=alt.Bin(maxbins=50), title="FDR-adjusted p-value"),
        y=alt.Y("count():Q", title="Count")
    ).properties(
        width=500,
        height=200,
        title={
            "text": "Distribution of FDR-adjusted p-values",
            "subtitle": f"Significant hits at 1% FDR highlighted (n={n_sig})",
            "anchor": "start",
        }
    )

    (pval_chart & fdr_chart).resolve_scale(y='independent')
    return


@app.cell
def _(df_sig):
    df_sig
    return


@app.cell
def _(alt, results_flt):
    df_sig = results_flt[results_flt["significant"]]

    df_count = (
        df_sig.groupby(["response_id", "day"])
           .size()
           .reset_index(name="n_hits")
    )

    chart = (
        alt.Chart(df_count)
        .mark_bar()
        .encode(
            x=alt.X("response_id:N", sort='-y', title="gRNA Target"),
            xOffset="day:N",  # <-- key for grouped bar charts
            y=alt.Y("n_hits:Q", title="# Significant Interactions"),
            color=alt.Color("day:N", title="Day"),
            tooltip=["response_id", "day", "n_hits"]
        )
        .properties(
            width=600,
            height=350,
            title="Significant interactions per gRNA target across days"
        )
        .configure_axis(labelFontSize=12, titleFontSize=14)
        .configure_title(fontSize=16, anchor="start")
    )

    chart
    return (df_sig,)


@app.cell
def _(results_flt):
    results_flt.sample(1000)
    return


if __name__ == "__main__":
    app.run()
