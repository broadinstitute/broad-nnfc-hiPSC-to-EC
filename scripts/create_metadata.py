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
    return (pd,)


@app.cell
def _(mo):
    mo.md(r"""### Create metadata""")
    return


@app.cell
def _(pd):
    # Add channel info (join instead of per-row .loc)
    channel_df = pd.read_csv("metadata/hiPSC-EC/channel_info.csv")

    # Extract metadata from Sample IDs and channel_info.csv
    pattern = r'^R(?P<replicate>\d+)-D(?P<day>\d+)-(?P<chip>[A-Za-z]+)-[Bb](?P<batch>\d+)$'
    parts = channel_df['sample'].str.extract(pattern)

    # Validate sample name structure
    bad = channel_df.loc[parts.isna().any(axis=1), 'sample']
    if len(bad):
        raise ValueError(f"Sample IDs not matching R#-D#-<Chip>-b#: {bad.tolist()}")

    # types
    parts = parts.astype({'replicate': 'int64', 'day': 'int64', 'batch': 'int64'})

    metadata_df = (
        channel_df[['sample']]
        .join(parts)
        .merge(channel_df, on='sample', how='left')
        .rename(columns={'channel':'channel'})
    )
    # warn if any channel missing
    missing_ch = metadata_df['channel'].isna()
    if missing_ch.any():
        print("WARNING: missing Channel for Samples:", metadata_df.loc[missing_ch, 'sample'].tolist())

    # 3) labels (vectorized; adjust zero padding to taste)
    metadata_df['sample_label_short'] = (
        'Rep_' + metadata_df['replicate'].astype(str) +
        '_Day_' + metadata_df['day'].astype(str)
    )

    metadata_df['sample_label_long'] = (
        'Replicate ' + metadata_df['replicate'].astype(str) +
        ' Day ' + metadata_df['day'].astype(str) +
        ' Chip ' + metadata_df['chip'] +
        ' Batch ' + metadata_df['batch'].astype(str) +
        ' Channel ' + metadata_df['channel'].astype(str)
    )

    # optional: a machine-safe id
    metadata_df['sample_id'] = (
        'rep' + metadata_df['replicate'].astype(str).str.zfill(2) + '_' +
        'day'   + metadata_df['day'].astype(str).str.zfill(2) + '_' +
        'chip_' + metadata_df['chip'].str.lower() + '_' +
        'batch' + metadata_df['batch'].astype(str).str.zfill(2) + '_' +
        'channel' + metadata_df['channel'].astype(str).str.zfill(2)
    )

    # save metadata
    metadata_df.to_csv('metadata/hiPSC-EC/metadata_sample.csv', index=False)
    metadata_df
    return


if __name__ == "__main__":
    app.run()
