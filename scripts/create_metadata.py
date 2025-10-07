import marimo

__generated_with = "0.15.2"
app = marimo.App()


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell
def _():
    import io
    import json
    import re
    import altair as alt
    import pandas as pd
    from functools import partial
    from pathlib import Path
    from google.cloud import storage
    return io, json, partial, pd, re, storage


@app.cell
def _(mo):
    mo.md(r"""### Create metadata""")
    return


@app.cell(hide_code=True)
def _(io, pd, storage):
    def read_gcs_file_to_pandas(bucket_name: str, file_path: str) -> pd.DataFrame:
        """
        Reads a file from GCS into a pandas DataFrame using an in-memory buffer.

        Args:
            bucket_name (str): Name of the GCS bucket.
            file_path (str): Full path to the file (blob name).

        Returns:
            pandas.DataFrame: The DataFrame created from the GCS file.
        """
        # 1. Initialize the GCS Client
        storage_client = storage.Client()

        # 2. Get the Blob (file) object
        bucket = storage_client.bucket(bucket_name)
        blob = bucket.blob(file_path)

        # 3. Download the content to an in-memory buffer (bytes)
        # The 'download_as_bytes()' method is efficient for this.
        file_bytes = blob.download_as_bytes()

        # 4. Create a BytesIO object (the file handler)
        # Pandas can read directly from this file-like object
        file_handler = io.BytesIO(file_bytes)

        # 5. Read the buffer into a DataFrame
        # Replace pd.read_csv with pd.read_parquet if needed
        df = pd.read_csv(file_handler)

        return df
    return (read_gcs_file_to_pandas,)


@app.cell(hide_code=True)
def _(partial, pd, read_gcs_file_to_pandas):
    # Load the FASTQ list
    bucket_name = "250806_d0_250617_ec-large-scale"
    fastq_paths = [
        "WALKUP-19009/250923_SL-EXL_0267_B232FYWLT4/2486_1758998285/fastq/Reports/fastq_list.csv",
        "WALKUP-18819/250821_SL-EXP_0163_A232FMMLT4/2311_1755887025/fastq/Reports/fastq_list.csv"
    ]

    # Use functools.partial to fix the bucket_name argument
    # This creates a new function that only needs the file_path
    read_with_bucket = partial(read_gcs_file_to_pandas, bucket_name)

    # Use a list comprehension to execute the function for all paths
    # The result is a list of pandas DataFrames
    list_of_dfs = [read_with_bucket(path) for path in fastq_paths]

    # Concatenate all DataFrames into a single DataFrame
    fastq_df = pd.concat(list_of_dfs, ignore_index=True)

    print("Resulting object type:", type(fastq_df))
    print(f"Total rows in final DataFrame: {len(fastq_df)}")

    fastq_df
    return (fastq_df,)


@app.cell
def _(fastq_df, re):
    # Experiment independent clean up

    # Extract the short sample name
    fastq_df['sample'] = fastq_df['RGSM'].str.split('_').str[-1].str[:-2]

    # Add an experiment name based on the name of the files :(
    fastq_df['experiment'] = fastq_df.apply(
        lambda row: 'wtc11' if row['sample'].startswith('1') else "endothelial_differentiation",
        axis=1
    )

    # Assign library_type to gene_expression or crispr_guide based if RGSM strings ends in m or g
    fastq_df['library_type'] = fastq_df['RGSM'].apply(
        lambda x: 'gene_expression' if x.endswith('m') else "crispr_guide"
    )

    # Update the FASTQ paths
    # Create the mapping
    path_map = {
        "/seq/dragen/fastqs/prod/SL-EXP/250821_SL-EXP_0163_A232FMMLT4/dragen/2025-08-21--15-00-17/fastq/":
            "gs://250806_d0_250617_ec-large-scale/WALKUP-18819/250821_SL-EXP_0163_A232FMMLT4/2311_1755887025/fastq/",
        "seq/dragen/fastqs/prod/SL-EXL/250923_SL-EXL_0267_B232FYWLT4/dragen/2025-09-27--08-55-54/fastq/":
            "gs://250806_d0_250617_ec-large-scale/WALKUP-19009/250923_SL-EXL_0267_B232FYWLT4/2486_1758998285/fastq/",
    }

    file_columns = ['Read1File', 'Read2File']

    # Build a prefix regex map:
    prefix_map = {
        rf'^/?{re.escape(src.lstrip("/"))}': dst
        for src, dst in path_map.items()
    }

    # Do substring replacement (vectorized) on just those columns
    fastq_df.loc[:, file_columns] = fastq_df.loc[:, file_columns].replace(prefix_map, regex=True)

    fastq_df
    return


@app.cell
def _(fastq_df):
    # WTC11 sample metadata creation
    metadata_wtc11 = fastq_df[fastq_df['experiment'] == "wtc11"].copy()

    # These are fixed because we only have one replicate and one batch
    metadata_wtc11["batch"] = 1
    metadata_wtc11["replicate"] = 1

    # Extract the metadata from the sample name :(
    # Create regex
    pattern_wtc11 = r'^(?P<viral_volume>\d+)-D(?P<day>\d+)-r(?P<round>\d+)-(?P<chip>[A-Z])(?P<channel>\d+)$'

    # Apply regex
    wtc11_sample_parse = metadata_wtc11['sample'].str.extract(pattern_wtc11)
    wtc11_sample_parse = wtc11_sample_parse[["day", "chip", "channel", "viral_volume", "round"]]
    # Merge results back to the dataframe
    metadata_wtc11 = (
        metadata_wtc11
        .join(wtc11_sample_parse)
    )

    metadata_wtc11.reset_index


    metadata_wtc11['sample_label_long'] = (
        'Replicate ' + metadata_wtc11['replicate'].astype(str) +
        ' Day ' + metadata_wtc11['day'].astype(str) +
        ' Chip ' + metadata_wtc11['chip'] +
        ' Batch ' + metadata_wtc11['batch'].astype(str) +
        ' Channel ' + metadata_wtc11['channel'].astype(str) + 
        ' Viral Volume ' + metadata_wtc11['viral_volume'].astype(str)
    )

    # optional: a machine-safe id
    metadata_wtc11['sample_id'] = (
        'rep' + metadata_wtc11['replicate'].astype(str).str.zfill(2) + '_' +
        'day'   + metadata_wtc11['day'].astype(str).str.zfill(2) + '_' +
        'chip_' + metadata_wtc11['chip'].str.lower() + '_' +
        'batch' + metadata_wtc11['batch'].astype(str).str.zfill(2) + '_' +
        'channel' + metadata_wtc11['channel'].astype(str).str.zfill(2) + '_' +
        'viral_volume' + metadata_wtc11['viral_volume'].astype(str)
    )

    metadata_wtc11 = metadata_wtc11[
        ['sample_id', 'sample_label_long'] +  metadata_wtc11.columns.tolist()[:-2]
        ]

    metadata_wtc11.to_csv("metadata/hiPSC-HPC/sample_metadata_wtc11.tsv", sep="\t", index=False)

    metadata_wtc11
    return (metadata_wtc11,)


@app.cell
def _(fastq_df, pd):
    # Endothelial differentiation metadata creation
    metadata_endothelial = fastq_df[fastq_df['experiment'] == "endothelial_differentiation"].copy()

    # Add channel info because they are missing even in the sample name
    channel_df = pd.read_csv("metadata/hiPSC-EC/channel_info.csv")

    # Extract the metadata from the sample name :(
    # Create regex
    pattern_endothelial = r'^R(?P<replicate>\d+)-D(?P<day>\d+)-(?P<chip>[A-Za-z]+)-[Bb](?P<batch>\d+).*$'
    # Apply regex
    endothelial_sample_parse = metadata_endothelial['sample'].str.extract(pattern_endothelial)
    endothelial_sample_parse = endothelial_sample_parse[["batch", "replicate", "day", "chip"]]
    # Merge results back to the dataframe
    metadata_endothelial = (
        metadata_endothelial
        .join(endothelial_sample_parse)
        .merge(channel_df, on='sample', how='left')
        .rename(columns={'channel':'channel'})
    )

    metadata_endothelial.reset_index

    metadata_endothelial['sample_label_long'] = (
        'Replicate ' + metadata_endothelial['replicate'].astype(str) +
        ' Day ' + metadata_endothelial['day'].astype(str) +
        ' Chip ' + metadata_endothelial['chip'] +
        ' Batch ' + metadata_endothelial['batch'].astype(str) +
        ' Channel ' + metadata_endothelial['channel'].astype(str)
    )

    # optional: a machine-safe id
    metadata_endothelial['sample_id'] = (
        'rep' + metadata_endothelial['replicate'].astype(str).str.zfill(2) + '_' +
        'day'   + metadata_endothelial['day'].astype(str).str.zfill(2) + '_' +
        'chip_' + metadata_endothelial['chip'].str.lower() + '_' +
        'batch' + metadata_endothelial['batch'].astype(str).str.zfill(2) + '_' +
        'channel' + metadata_endothelial['channel'].astype(str).str.zfill(2)
    )

    metadata_endothelial = metadata_endothelial[
        ['sample_id', 'sample_label_long'] +  metadata_endothelial.columns.tolist()[:-2]
        ]

    metadata_endothelial.to_csv("metadata/hiPSC-EC/sample_metadata_endothelial.tsv", sep="\t", index=False)

    metadata_endothelial
    return (metadata_endothelial,)


@app.cell
def _(json, pd, re):
    def build_terra_fastq_table(
        df: pd.DataFrame,
        *,
        sample_col: str = "sample_id",
        library_col: str = "library_type",
        r1_col: str = "Read1File",
        r2_col: str = "Read2File",
        # map library_type -> (out_R1_col, out_R2_col)
        library_map: dict = None,  # defaults below
        exclude_cols: set = {"Lane", "RGSM", "RGID"},
        lane_regex: str = r"_L(\d{3})_",
        rename_for_terra: bool = True,
        out_tsv_path: str | None = None,
        warn_on_inconsistent: bool = True,
    ) -> pd.DataFrame:
        """
        Build a Terra-friendly per-sample table:
          - keeps all non-FASTQ columns identical within sample (first value)
          - creates JSON-array string columns per library type for R1/R2
          - drops columns in `exclude_cols` before merging

        Returns the resulting DataFrame and optionally writes a TSV.
        """
        if library_map is None:
            library_map = {
                "gene_expression": ("cdna_read1_files", "cdna_read2_files"),
                "crispr_guide":   ("crispr_read1_files", "crispr_read2_files"),
            }

        # --- working copy (drop excluded columns) ---
        df = df.drop(columns=[c for c in exclude_cols if c in df.columns]).copy()

        # --- helpers (sorting & aggregation) ---
        LANE_RX = re.compile(lane_regex)
        def lane_number(u: str) -> int:
            m = LANE_RX.search(u or "")
            return int(m.group(1)) if m else 10**6

        def collect_sorted_unique(series: pd.Series) -> list[str]:
            vals = [str(x) for x in series.dropna()]
            seen, out = set(), []
            for v in vals:
                if v not in seen:
                    seen.add(v); out.append(v)
            out.sort(key=lambda u: (lane_number(u), u))
            return out

        # --- figure out which columns to carry through unchanged ---
        group_cols = [sample_col]
        fastq_cols = [r1_col, r2_col, library_col]
        other_cols = [c for c in df.columns if c not in set(group_cols + fastq_cols)]

        # --- optional validation: identical per sample ---
        if warn_on_inconsistent and other_cols:
            nunq = df.groupby(sample_col)[other_cols].nunique(dropna=False)
            bad_samples = nunq.gt(1).any(axis=1)
            if bad_samples.any():
                print(
                    "WARNING: some non-FASTQ columns vary within a sample:",
                    bad_samples.index[bad_samples].tolist()
                )

        # --- aggregate FASTQs per (sample, library) ---
        agg = (
            df.groupby([sample_col, library_col], as_index=False, sort=False)
              .agg(**{
                  r1_col: (r1_col, collect_sorted_unique),
                  r2_col: (r2_col, collect_sorted_unique),
              })
        )

        # --- pivot wide: MultiIndex columns (r1/r2 x library_type) ---
        wide = agg.pivot(index=sample_col, columns=library_col, values=[r1_col, r2_col])

        def safe_list(key):
            s = wide.get(key)
            if s is None:
                return pd.Series([[] for _ in range(len(wide))], index=wide.index)
            return s.apply(lambda x: x if isinstance(x, list) else ([] if pd.isna(x) else [str(x)]))

        # --- build output FASTQ columns per library_map ---
        fastq_wide = pd.DataFrame(index=wide.index)
        for lib, (out_r1, out_r2) in library_map.items():
            fastq_wide[out_r1] = safe_list((r1_col, lib))
            fastq_wide[out_r2] = safe_list((r2_col, lib))
        fastq_wide = fastq_wide.reset_index()

        # --- keep other metadata unchanged (first per sample) ---
        if other_cols:
            meta_first = df.groupby(sample_col, as_index=False)[other_cols].first()
        else:
            meta_first = df[[sample_col]].drop_duplicates()

        # --- combine & JSON-encode arrays for Terra ---
        out = meta_first.merge(fastq_wide, on=sample_col, how="left")

        json_cols = [c for c in out.columns if c not in other_cols + [sample_col]]
        for col in json_cols:
            out[col] = out[col].apply(lambda lst: json.dumps(lst or [], ensure_ascii=False))

        # --- rename entity id for Terra ---
        if rename_for_terra:
            out = out.rename(columns={sample_col: f"entity:{sample_col}"})

        # --- optional write ---
        if out_tsv_path:
            out.to_csv(out_tsv_path, sep="\t", index=False)

        return out
    return (build_terra_fastq_table,)


@app.cell
def _(build_terra_fastq_table, metadata_endothelial, metadata_wtc11):
    # For your current df
    terra_df1 = build_terra_fastq_table(
        metadata_wtc11,
        out_tsv_path="metadata/hiPSC-HPC/sample_metadata_wtc11_terra.tsv"
    )
    terra_df2 = build_terra_fastq_table(
        metadata_endothelial,
        out_tsv_path="metadata/hiPSC-EC/sample_metadata_endothelial_terra.tsv"
    )
    return


if __name__ == "__main__":
    app.run()
