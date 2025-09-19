import marimo

__generated_with = "0.15.2"
app = marimo.App()


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell
def _():
    import altair as alt
    import pandas as pd
    return alt, pd


@app.cell
def _(mo):
    mo.md(r"""# Initial QC for the hiPSC to EC samples""")
    return


@app.cell
def _(mo):
    mo.md(r"""### Load metadata""")
    return


@app.cell(hide_code=True)
def _(pd):
    metadata_df = pd.read_csv('metadata/hiPSC-EC/metadata.csv')
    metadata_df
    return (metadata_df,)


@app.cell
def _(mo):
    mo.md(r"""### Cellranger metrics""")
    return


@app.cell(hide_code=True)
def _(metadata_df, pd):
    metrics_df = pd.read_table('data/ec-differentiation/cellranger_qc/metrics_qc_custom_mqc.tsv', sep='\t', comment="#")

    metrics_df.rename(columns={'Sample': 'sample'}, inplace=True)
    metrics_df = metrics_df.merge(
        metadata_df[['sample', 'sample_label_short', 'batch','chip', 'channel', 'replicate', 'day']],
        on='sample',
        how='left'
    )
    # Sanitize columns
    metrics_df.columns = metrics_df.columns.str.replace(':', '_', regex=False)
    metrics_df.columns = metrics_df.columns.str.replace(' ', '_', regex=False)

    metrics_df
    return (metrics_df,)


@app.cell(hide_code=True)
def _(alt, metrics_df):
    _chart = (
        alt.Chart(metrics_df)
        .mark_bar()
        .encode(
            x=alt.X(field='day', type='nominal'),
            y=alt.Y(field='replicate', type='nominal', aggregate='count'),
            color=alt.Color(field='replicate', type='nominal'),
            tooltip=[
                alt.Tooltip(field='day', format=',.0f', title=''),
                alt.Tooltip(field='replicate',
                            aggregate='count', 
                            format=',.0f', title=''),
                alt.Tooltip(field='replicate')
            ]
        )
        .properties(
            title='',
            height=290,
            width='container',
            config={
                'axis': {
                    'grid': False
                }
            }
        )
    )
    _chart
    return


@app.cell(hide_code=True)
def _(alt, metadata_df, mo):
    # Using the metadata_df plot the sample name using 
    # Chip as rows and Channel as column and facet by Batch
    palette= ["#f58518","#B7D7A8","#4c78a8","#e45756"]

    _chart = alt.Chart(metadata_df).mark_text(fontWeight='bold').encode(
        x='channel:N',
        y='chip:N',
        text='sample_label_short:N',
        color=alt.Color('replicate:N',
                        scale=alt.Scale(range=palette),
                        legend = None
                       )
    ).properties(
        width=500,
        height=300
    ).facet(
        row='batch:N'
    ).resolve_scale(
        x='independent'
    ).configure_axis(
        labelFontSize=12,
        titleFontSize=14,
        labelAngle=0
    ).configure_header(
        labelFontSize=14,
        titleFontSize=16
    )
    mo.ui.altair_chart(_chart)
    return


@app.cell(hide_code=True)
def _(metadata_df, metrics_df, mo):
    # Marimo ui element to select the metric to plot
    metric_selector = mo.ui.dropdown(options=metrics_df.columns.tolist(),
                                     label="Choose metric",
                                     value=metrics_df.columns.tolist()[1]
                                    )
    # Marimo ui multiselect to choose the replicate
    options_replicate = metadata_df["replicate"].unique().tolist()
    multiselect_replicate = mo.ui.multiselect(
        options=options_replicate,
        label="Replicate",
        value=options_replicate,
    )
    # Marimo ui multiselect to choose the day
    options_day = metadata_df["day"].unique().tolist()
    multiselect_day = mo.ui.multiselect(
        options=options_day,
        label="Day",
        value=options_day,
    )
    return (metric_selector,)


@app.cell(hide_code=True)
def _(alt, metric_selector, metrics_df, mo):
    # Formats available
    #',d' → 23,000,000
    #',.1f' → 23,000,000.0
    #'.1%' → 89.7% (if your data are fractions 0–1)
    #'~s' → 23M (SI shorthand)

    _chart = alt.Chart(metrics_df).mark_rect(stroke="white").encode(
        x='channel:N',
        y='chip:N',
        color=alt.Color(f"{metric_selector.value}:Q",
                        scale=alt.Scale(scheme='viridis')
                       ),
        tooltip=[
            alt.Tooltip('sample_label_short:N', title='Sample'),
            alt.Tooltip(
                field=metric_selector.value,
                type='quantitative',
                format=",d",
                title=metric_selector.value
            )
        ]
    ).properties(
        width=450,
        height=280
    ).facet(
        row='batch:N'
    ).resolve_scale(
        x='independent'
    ).configure_axis(
        labelFontSize=12,
        titleFontSize=14,
        labelAngle=0
    ).configure_header(
        labelFontSize=14,
        titleFontSize=16
    )
    mo.vstack([metric_selector, _chart])
    return


@app.cell(hide_code=True)
def _(alt, metric_selector, metrics_df, mo):
    # Combine the day and replicate selector together
    _brush = alt.selection_point(encodings=["x","color","shape"])

    # Plot histogram of the metric selected with the metric_selector
    # Faceted by batch, colored by replicate, and line type based on day
    _chart = (
            alt.Chart(metrics_df)
            .mark_point(
                filled=True,
                size=30
            ).encode(
                x=alt.X(
                    'day:N',
                    title='Day',
                    axis=alt.Axis(labelAngle=0),
                ),
                xOffset=alt.XOffset('replicate:N'),
                y=alt.Y(
                    f'{metric_selector.value}:Q',
                    title=metric_selector.value,
                    axis=alt.Axis(format=",d"),
                    scale=alt.Scale(zero=False),
                ),
                color=alt.Color('replicate:N', title='replicate'),
                shape=alt.Shape('batch:N', title="Batch"),
            ).properties(
        width=450,
        height=250
    ).add_params(_brush)
    )
    chart = mo.ui.altair_chart(
        _chart,
    )
    return (chart,)


@app.cell(hide_code=True)
def _(chart, mo):
    mo.vstack([chart, mo.ui.table(chart.value)])

    return


@app.cell
def _(metrics_df):
    metrics_df.groupby(['day','replicate'])['Estimated_Number_of_Cells'].sum()
    return


if __name__ == "__main__":
    app.run()
