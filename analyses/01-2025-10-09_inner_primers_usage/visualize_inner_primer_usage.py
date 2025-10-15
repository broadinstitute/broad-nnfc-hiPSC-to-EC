import marimo

__generated_with = "0.15.2"
app = marimo.App()


@app.cell
def _():
    import marimo as mo
    return


@app.cell
def _():
    import pandas as pd
    import numpy as np
    from collections import Counter, defaultdict
    import altair as alt
    return alt, pd


@app.cell
def _():
    output_folder = "results/01-2025-10-09_inner_primers_usage"
    inner_primer_usage = "gs://fc-76565551-28c3-4b6d-a048-e272103bcbd1/submissions/intermediates/a5f262a4-2e67-4187-8afa-b4ba48adc7e8/CheckInnerPrimerUsage/4661facf-682a-49f9-87c1-41e051a4e715/call-AnalyzePrimerUsage/rep01_day00_chip_a_batch01_channel01_detailed.tsv"
    summary_report = "gs://fc-76565551-28c3-4b6d-a048-e272103bcbd1/submissions/intermediates/a5f262a4-2e67-4187-8afa-b4ba48adc7e8/CheckInnerPrimerUsage/4661facf-682a-49f9-87c1-41e051a4e715/call-AnalyzePrimerUsage/rep01_day00_chip_a_batch01_channel01.txt"
    return (inner_primer_usage,)


@app.cell
def _(inner_primer_usage, pd):
    df2 = pd.read_csv(inner_primer_usage, sep="\t")
    df2
    return (df2,)


@app.cell
def _(df2):
    eff = (
        df2.assign(is_correct=df2['Match_Status'].eq('correct'))
           .groupby('Inner_Primer_Found')['is_correct']
           .mean()
           .rename('frac_correct')          # in [0,1]
    )
    # If you want percent:
    eff_pct = (eff * 100).rename('pct_correct')
    return (eff_pct,)


@app.cell
def _(eff_pct):
    eff_pct
    return


@app.cell
def _():
    return


@app.cell
def _(alt, eff_pct):
    chart =(
        alt.Chart(eff_pct)
            .mark_bar()
            .encode(
                x=alt.X("pct_correct:Q", bin=alt.Bin(maxbins=30), title="Primer efficiency (%)"),
                y=alt.Y("count():Q", title="# primers"),
                tooltip=[alt.Tooltip("count():Q", title="# primers")]
            )
            .properties(width=600, height=320, title="Distribution of primer efficiencies")
    )

    chart
    return


@app.cell
def _(df, pd):
    # 1) Map each primer to its intended gene (expect 1:1)
    #    If a primer maps to >1 intended gene in the data, we keep the most frequent.
    intended_map = (
        df.groupby(['Inner_Primer_Found', 'Intended_Target_Gene_ID'])
          .size().reset_index(name='n')
          .sort_values(['Inner_Primer_Found','n'], ascending=[True, False])
          .drop_duplicates('Inner_Primer_Found')
          .set_index('Inner_Primer_Found')['Intended_Target_Gene_ID']
          .to_dict()
    )

    # 2) Build the counts matrix (primer x observed gene)
    pivot = pd.crosstab(df['Inner_Primer_Found'], df['Observed_Gene_ID'])

    # 3) Order rows/cols so the expected gene for each primer lies on the diagonal
    row_order = sorted(pivot.index.tolist())  # or order by total counts if you prefer
    expected_cols_in_order = [intended_map.get(p, None) for p in row_order]
    expected_cols_in_order = [c for c in expected_cols_in_order if c is not None]

    other_cols = [c for c in pivot.columns if c not in set(expected_cols_in_order)]
    col_order = expected_cols_in_order + sorted(other_cols)

    pivot = pivot.reindex(index=row_order, columns=col_order, fill_value=0)
    return


if __name__ == "__main__":
    app.run()
