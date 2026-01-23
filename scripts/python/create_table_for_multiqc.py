import marimo

__generated_with = "0.18.4"
app = marimo.App()


@app.cell
def _():
    import marimo as mo
    return


@app.cell
def _():
    import pandas as pd
    import glob
    import os
    return glob, os, pd


@app.cell
def _(pd):
    pd.read_csv("../../results/2025-12-18_ec_multiqc_batch4/metrics2_qc/rep03_day00_chip_a_batch04_channel01_metrics_summary.csv").head()
    return


@app.cell
def _():
    out_path = "../../results/2025-12-18_ec_multiqc_batch4/metrics_qc_custom_mqc.tsv"
    return (out_path,)


@app.cell
def _(glob, os, out_path, pd):
    # --- Function to read, label, and clean the DataFrame ---

    def sanitize_and_label_dataframe(df):
        """Sanitizes all string columns by removing separators and percentages."""

        # Regex pattern to match commas (,), tabs (\t), spaces (\s), and percent signs (%)
        chars_to_remove = r'[,\t\s%]'

        # Select all columns with string data (dtype 'object')
        str_cols = df.select_dtypes(include=['object']).columns

        # Apply the replacement vectorially across all selected string columns
        for col in str_cols:
            # Use str.replace with regex=True for global replacement
            df[col] = df[col].str.replace(chars_to_remove, '', regex=True)

        return df

    # --- 1. Ingestion and Labeling ---

    def read_and_label_csv(file_path):
        """Reads a CSV, adds a 'source_file' column, and cleans the data."""
        # Read the file
        df = pd.read_csv(file_path)

        # Add the new column with the filename
        df['source_file'] = os.path.basename(file_path)

        # Sanitize the DataFrame before concatenating
        return sanitize_and_label_dataframe(df)


    # Ingestion loop: Read and process all files
    file_paths = glob.glob("../../results/2025-12-18_ec_multiqc_batch4/metrics2_qc/*_metrics_summary.csv")

    list_of_dfs = [read_and_label_csv(f) for f in file_paths]

    # Concatenate all DataFrames
    df = pd.concat(list_of_dfs, ignore_index=True)

    # --- 2. Final Export ---

    # Export the clean DataFrame to the output path
    df.to_csv(out_path, sep="\t", index=False)
    return


@app.cell
def _(out_path, pd):
    metrics = pd.read_csv(out_path, sep="\t")

    metrics['clean_name'] = metrics["source_file"].str.replace(
        "_metrics_summary.csv", "", regex=False
    )

    # 2. Define the extraction pattern to pull out named components
    # This pattern extracts the digits/letters following the descriptive prefix (rep, day, chip, batch, channel)
    extraction_pattern = (
        r'rep(?P<replicate>\d+)_day(?P<day>\d+)_chip_'
        r'(?P<chip_letter>[a-z]+)_batch(?P<batch_id>\d+)_channel(?P<channel_id>\d+)'
    )

    # 3. Apply the extraction and assign the results to new columns
    metrics[[
        'replicate',      # e.g., '02'
        'day',            # e.g., '02'
        'chip_letter',    # e.g., 'a'
        'batch_id',       # e.g., '03'
        'channel_id'      # e.g., '06'
    ]] = metrics['clean_name'].str.extract(extraction_pattern)

    # You can now drop the temporary 'clean_name' column if desired
    metrics = metrics.drop(columns=['clean_name'])
    metrics.head()
    return (metrics,)


@app.cell
def _(metrics):
    # Compute cells recovered per day with number of cells and median reads per cell
    metrics.groupby(['day'])[['Estimated Number of Cells', 'Mean Reads per Cell']].sum()
    return


@app.cell
def _(metrics):
    metrics.groupby('day').agg({
        'Estimated Number of Cells': 'sum', 
        'Mean Reads per Cell': 'mean'
    })
    return


if __name__ == "__main__":
    app.run()
