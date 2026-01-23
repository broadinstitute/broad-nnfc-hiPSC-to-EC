import marimo

__generated_with = "0.17.7"
app = marimo.App()


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell
def _():
    import re
    import pandas as pd
    import numpy as np
    from collections import Counter, defaultdict
    import altair as alt

    alt.data_transformers.disable_max_rows()
    return alt, pd, re


@app.cell
def _(mo):
    mo.md(r"""
    # Analysis of inner primer usage
    This notebook contains the results of the analysis of inner primer usage for the hiPSC -> EC and WTC11.
    In order for the results to be comparable, this notebook will only consider the `Day 0` timepoint for the hiPSC -> EC experiment.

    WARNING: All the statistics have been generated counting **reads** from the cellranger BAM alignment file. No UMI-based deduplication filter has been applied.
    """)
    return


@app.cell
def _():
    output_folder = "results/01-2025-10-09_inner_primers_usage"
    inner_primer_usage_ec_fnp = "gs://fc-76565551-28c3-4b6d-a048-e272103bcbd1/submissions/final-outputs/a5f262a4-2e67-4187-8afa-b4ba48adc7e8/CheckInnerPrimerUsage/4661facf-682a-49f9-87c1-41e051a4e715/call-AnalyzePrimerUsage/rep01_day00_chip_a_batch01_channel01_detailed.tsv"
    metadata_genes_ec_fnp = "../metadata/metadata_genes_final.csv"
    inner_primer_usage_wtc11_fnp = "gs://fc-76565551-28c3-4b6d-a048-e272103bcbd1/submissions/intermediates/50869f0a-cc33-4c65-bd42-9d7db75645d6/CheckInnerPrimerUsage/e7ae4018-ff29-46ac-994a-b87bdb033829/call-AnalyzePrimerUsage/rep01_day00_chip_a_batch01_channel01_volume100_detailed.tsv"
    metadata_genes_wtc11_fnp = "metadata/metadata_genes_with_ids.csv"
    return (
        inner_primer_usage_ec_fnp,
        inner_primer_usage_wtc11_fnp,
        metadata_genes_wtc11_fnp,
    )


@app.cell
def _(mo):
    mo.md(r"""
    ### Load results and metadata for the hiPSC -> EC
    """)
    return


@app.cell
def _(inner_primer_usage_ec_fnp, pd):
    inner_primer_usage_ec = pd.read_csv(inner_primer_usage_ec_fnp, sep="\t")
    inner_primer_usage_ec
    return (inner_primer_usage_ec,)


@app.cell
def _(pd):
    metadata_genes_ec = pd.read_csv("../../metadata/metadata_genes_final.csv")
    metadata_genes_ec
    return (metadata_genes_ec,)


@app.cell
def _(inner_primer_usage_ec):
    eff_ec = (
        inner_primer_usage_ec.assign(is_correct=inner_primer_usage_ec['Match_Status'].eq('correct'))
           .groupby('Inner_Primer_Found')['is_correct']
           .mean()
           .rename('frac_correct')          # in [0,1]
    )

    eff_ec_pct = (eff_ec * 100).rename('pct_correct')
    return (eff_ec_pct,)


@app.cell
def _(eff_ec_pct):
    # Save results to file
    eff_ec_pct.to_csv(f"../../results/01-2025-10-09_inner_primers_usage/inner_primer_efficiency_ec_pct.tsv", sep="\t", header=True)
    eff_ec_pct.head()
    return


@app.cell
def _(alt, eff_ec_pct, pd):
    # If eff_ec_pct is a list/Series/np.array, wrap it into a DataFrame:
    df = pd.DataFrame({"pct_correct": eff_ec_pct})

    chart1 = (
        alt.Chart(df)
        .mark_bar()
        .encode(
            x=alt.X("pct_correct:Q", bin=True, title="Primer Efficiency (%)"),
            y=alt.Y("count()", title="Frequency")
        )
        .properties(title="Histogram of Primer Efficiency")
    )

    chart1
    return


@app.cell
def _(mo):
    mo.md(r"""
    ### Load results and metadata for the WTC11
    """)
    return


@app.cell
def _(inner_primer_usage_wtc11_fnp, pd):
    inner_primer_usage_wtc11 = pd.read_csv(inner_primer_usage_wtc11_fnp, sep="\t")
    inner_primer_usage_wtc11
    return (inner_primer_usage_wtc11,)


@app.cell
def _(metadata_genes_wtc11_fnp, pd):
    metadata_genes_wtc11 = pd.read_csv(metadata_genes_wtc11_fnp)
    metadata_genes_wtc11
    return (metadata_genes_wtc11,)


@app.cell
def _(inner_primer_usage_wtc11):
    eff_wtc11 = (
        inner_primer_usage_wtc11.assign(is_correct=inner_primer_usage_wtc11['Match_Status'].eq('correct'))
           .groupby('Inner_Primer_Found')['is_correct']
           .mean()
           .rename('frac_correct')          # in [0,1]
    )

    eff_wtc11_pct = (eff_wtc11 * 100).rename('pct_correct')
    return eff_wtc11, eff_wtc11_pct


@app.cell
def _(eff_wtc11):
    eff_wtc11
    return


@app.cell
def _(mo):
    mo.md(r"""
    ### Compare efficiency for the shared dial-out genes
    The efficiency is defined as the fraction of UMIs with correct inner primer match out of all reads with that inner primer found.
    """)
    return


@app.cell
def _(metadata_genes_ec, metadata_genes_wtc11, mo):
    # Number of genes in common
    common_genes = set(metadata_genes_ec['intended_target_gene_id']).intersection(set(metadata_genes_wtc11['intended_target_gene_id']))
    mo.vstack([mo.md("Dial out genes in common:"),len(common_genes)])
    return (common_genes,)


@app.cell
def _(
    common_genes,
    eff_ec_pct,
    eff_wtc11_pct,
    metadata_genes_ec,
    metadata_genes_wtc11,
    pd,
):
    # Check if we used the same inner primers for the genes in common
    # Get output dataframe intended_target_gene_id, intended_target_gene_symbol inner_primer_ec, inner_primer_wtc11
    metadata_genes_ec_sub = metadata_genes_ec[metadata_genes_ec['intended_target_gene_id'].isin(common_genes)][['intended_target_gene_id', 'intended_target_gene_symbol', 'inner_primer']]
    metadata_genes_wtc11_sub = metadata_genes_wtc11[metadata_genes_wtc11['intended_target_gene_id'].isin(common_genes)][['intended_target_gene_id', 'intended_target_gene_symbol', 'inner_primer']]
    merged_metadata = pd.merge(metadata_genes_ec_sub, metadata_genes_wtc11_sub, on=['intended_target_gene_id', 'intended_target_gene_symbol'], suffixes=('_ec', '_wtc11'))
    merged_metadata['same_inner_primer'] = merged_metadata['inner_primer_ec'] == merged_metadata['inner_primer_wtc11']
    merged_metadata

    # Merge the efficiency data with the metadata
    eff_ec_df = eff_ec_pct.reset_index().rename(columns={'Inner_Primer_Found': 'inner_primer', 'pct_correct': 'pct_correct_ec'})
    eff_wtc11_df = eff_wtc11_pct.reset_index().rename(columns={'Inner_Primer_Found': 'inner_primer', 'pct_correct': 'pct_correct_wtc11'})
    merged_metadata = pd.merge(merged_metadata, eff_ec_df, left_on='inner_primer_ec', right_on='inner_primer', how='left')
    merged_metadata = pd.merge(merged_metadata, eff_wtc11_df, left_on='inner_primer_wtc11', right_on='inner_primer', how='left')
    merged_metadata = merged_metadata.drop(columns=['inner_primer_x', 'inner_primer_y'])
    merged_metadata
    return (merged_metadata,)


@app.cell
def _(merged_metadata, mo):
    mo.md(r"""We used the same inner primer for {} out of {} common dial-out genes.""".format(merged_metadata['same_inner_primer'].sum(), len(merged_metadata)))
    return


@app.cell
def _(merged_metadata):
    merged_metadata
    return


@app.cell
def _(alt, merged_metadata, mo):
    mo.md(r"""#### Scatter plot of primer efficiency""")
    chart = (
        alt.Chart(merged_metadata)
            .mark_circle(size=60)
            .encode(
                x=alt.X("pct_correct_ec:Q", title="Primer efficiency hiPSC->EC (%)"),
                y=alt.Y("pct_correct_wtc11:Q", title="Primer efficiency WTC11 (%)"),
                color=alt.Color("same_inner_primer:N", title="Same inner primer"),
                tooltip=[
                    alt.Tooltip("intended_target_gene_symbol:N", title="Gene symbol"),
                    alt.Tooltip("intended_target_gene_id:N", title="Gene ID"),
                    alt.Tooltip("inner_primer_ec:N", title="Inner primer hiPSC->EC"),
                    alt.Tooltip("pct_correct_ec:Q", title="Primer efficiency hiPSC->EC (%)"),
                    alt.Tooltip("inner_primer_wtc11:N", title="Inner primer WTC11"),
                    alt.Tooltip("pct_correct_wtc11:Q", title="Primer efficiency WTC11 (%)"),
                    alt.Tooltip("same_inner_primer:N", title="Same inner primer")
                ]
            )
            .properties(width=600, height=400, title="Primer efficiency comparison")
            .interactive()
    )
    chart
    return


@app.cell
def _(merged_metadata):
    merged_metadata
    return


@app.cell
def _(mo):
    mo.md(r"""
    If the genes that we don't expect are caused by a promiscous inner primer we can the information provided by the shared genes. For the 13 out of the 36 with the **shared inner_primer** we expect **the off target genes to be the same**. Conversely, for the 23 with different inner primers we expect different off target genes.
    """)
    return


@app.cell
def _(inner_primer_usage_ec):
    inner_primer_usage_ec.columns
    return


@app.cell
def _(inner_primer_usage_ec, inner_primer_usage_wtc11, merged_metadata, pd):
    # Take the inner primer for hiPSC -> EC that is shared with WTC11 and extract all the observed genes for it.
    shared_primers = merged_metadata[merged_metadata['same_inner_primer']]['inner_primer_ec'].unique().tolist()

    # Find the observed genes for these primers in the hiPSC -> EC data and WTC11
    temp_ec = (
        inner_primer_usage_ec
        .groupby(["Inner_Primer_Found", "Observed_Gene_ID"], dropna=False)
        .size()
        .reset_index(name="n")
        .rename(columns={"Inner_Primer_Found": "inner_primer", "Observed_Gene_ID": "observed_gene_id"})
    )
    observed_genes_ec_shared = temp_ec[
        temp_ec['inner_primer'].isin(shared_primers)
    ]
    observed_genes_ec_shared

    temp_wtc11 = (
        inner_primer_usage_wtc11
        .groupby(["Inner_Primer_Found", "Observed_Gene_ID"], dropna=False)
        .size()
        .reset_index(name="n")
        .rename(columns={"Inner_Primer_Found": "inner_primer", "Observed_Gene_ID": "observed_gene_id"})
    )
    observed_genes_wtc11_shared = temp_wtc11[
        temp_wtc11['inner_primer'].isin(shared_primers)
        ]
    # compute Jaccard index for each primer
    def jaccard_index(set1, set2):
        intersection = len(set1.intersection(set2))
        union = len(set1.union(set2))
        if union == 0:
            return 0.0
        return intersection / union

    jaccard_results = []
    for primer in shared_primers:
        genes_ec = set(observed_genes_ec_shared[observed_genes_ec_shared['inner_primer'] == primer]['observed_gene_id'].dropna().unique().tolist())
        genes_wtc11 = set(observed_genes_wtc11_shared[observed_genes_wtc11_shared['inner_primer'] == primer]['observed_gene_id'].dropna().unique().tolist())
        jaccard = jaccard_index(genes_ec, genes_wtc11)
        jaccard_results.append({
            'inner_primer': primer,
            'jaccard_index': jaccard,
            'n_genes_ec': len(genes_ec),
            'n_genes_wtc11': len(genes_wtc11),
            'n_shared_genes': len(genes_ec.intersection(genes_wtc11))
        })
    jaccard_results_df = pd.DataFrame(jaccard_results)
    jaccard_results_df
    return jaccard_index, jaccard_results_df


@app.cell
def _(mo):
    mo.md(r"""
    Now let's look at the 23 genes with differen inner primers
    """)
    return


@app.cell
def _(
    inner_primer_usage_ec,
    inner_primer_usage_wtc11,
    jaccard_index,
    merged_metadata,
    pd,
    re,
):
    # Choose only genes with different primers
    same_gene_different_primer = merged_metadata[~merged_metadata['same_inner_primer']]['intended_target_gene_id'].unique().tolist()
    # Strip gene version
    same_gene_different_primer = [
        re.sub(r"\.\d+$", "", gene_id) for gene_id in same_gene_different_primer
    ]
    temp2_ec = (
        inner_primer_usage_ec
        .groupby(["Intended_Target_Gene_ID", "Observed_Gene_ID"], dropna=False)
        .size()
        .reset_index(name="n")
        .rename(columns={"Intended_Target_Gene_ID": "intended_target_gene_id", "Observed_Gene_ID": "observed_gene_id"})
    )

    observed_genes_ec_unique = temp2_ec[
        temp2_ec['intended_target_gene_id'].isin(same_gene_different_primer)
    ]

    temp2_wtc11 = (
        inner_primer_usage_wtc11
        .groupby(["Intended_Target_Gene_ID", "Observed_Gene_ID"], dropna=False)
        .size()
        .reset_index(name="n")
        .rename(columns={"Intended_Target_Gene_ID": "intended_target_gene_id", "Observed_Gene_ID": "observed_gene_id"})
    )
    observed_genes_wtc11_unique = temp2_wtc11[
        temp2_wtc11['intended_target_gene_id'].isin(same_gene_different_primer)
        ]

    jaccard_results_different = []
    for shared_gene in same_gene_different_primer:
        genes_ec2 = set(observed_genes_ec_unique[observed_genes_ec_unique['intended_target_gene_id'] == shared_gene]['observed_gene_id'].dropna().unique().tolist())
        genes_wtc112 = set(observed_genes_wtc11_unique[observed_genes_wtc11_unique['intended_target_gene_id'] == shared_gene]['observed_gene_id'].dropna().unique().tolist())
        jaccard2 = jaccard_index(genes_ec2, genes_wtc112)
        jaccard_results_different.append({
            'intended_target_gene_id': shared_gene,
            'jaccard_index': jaccard2,
            'n_genes_ec': len(genes_ec2),
            'n_genes_wtc11': len(genes_wtc112),
            'n_shared_genes': len(genes_ec2.intersection(genes_wtc112))
        })
    jaccard_results_different_df = pd.DataFrame(jaccard_results_different)
    jaccard_results_different_df
    return (jaccard_results_different_df,)


@app.cell
def _(jaccard_results_different_df):
    jaccard_results_different_df.melt(id_vars=['intended_target_gene_id'], value_vars=['jaccard_index'], var_name='metric', value_name='value')
    return


@app.cell
def _(alt, jaccard_results_df, jaccard_results_different_df, pd):
    # Plot in the same plot two different columns the jitter of the jaccard results for jaccard_results_different_df and jaccard_results
    jaccard_results_df_melted = jaccard_results_df.melt(id_vars=['inner_primer'], value_vars=['jaccard_index'], var_name='metric', value_name='value')
    jaccard_results_df_melted
    jaccard_results_different_df_melted = jaccard_results_different_df.melt(id_vars=['intended_target_gene_id'], value_vars=['jaccard_index'], var_name='metric', value_name='value')
    jaccard_results_df_melted['type'] = 'Same inner primer'
    jaccard_results_different_df_melted['type'] = 'Different inner primer'
    jaccard_combined = pd.concat([jaccard_results_df_melted, jaccard_results_different_df_melted], ignore_index=True)
    jaccard_combined
    chart_jaccard = (
        alt.Chart(jaccard_combined)
        .mark_circle(size=60)
        .encode(
            x=alt.X("type:N", title=""),
            y=alt.Y("value:Q", title="Value"),
            color=alt.Color("type:N", title="Inner primer type"),
            tooltip=[
                alt.Tooltip("inner_primer:N", title="Inner primer"),
                alt.Tooltip("intended_target_gene_id:N", title="Gene ID"),
                alt.Tooltip("metric:N", title="Metric"),
                alt.Tooltip("value:Q", title="Value"),
                alt.Tooltip("type:N", title="Inner primer type")
            ]
        )
        .properties(width=600, height=400, title="Jaccard index comparison")
        .interactive()
    )
    chart_jaccard
    return


@app.cell
def _(mo):
    mo.md(r"""
    ### hiPSC -> EC
    Analysis specific to the hiPSC -> EC experiment.
    """)
    return


if __name__ == "__main__":
    app.run()
