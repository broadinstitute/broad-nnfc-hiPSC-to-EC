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


@app.cell(hide_code=True)
def _():
    import uuid
    import random


    def _get_accession_type_abbr(accession_type: str) -> str:
        """
        Returns the abbreviation for a given accession type.
        Raises ValueError if the accession type is invalid.
        """
        valid_types = {'guide': "GU", 'element': "EL", 'outer_primer': "OP", 'inner_primer': "IP"}
        if accession_type.lower() not in valid_types:
            raise ValueError(f"Invalid accession_type '{accession_type}'. Must be one of {valid_types}.")
        return valid_types[accession_type.lower()]


    def _get_custom_namespace(accession_type: str) -> uuid.UUID:
        """
        Returns a custom namespace UUID based on the accession type.
        """
        NAMESPACE_DNS = uuid.NAMESPACE_DNS
        namespace = uuid.uuid5(NAMESPACE_DNS, "hgrm." + accession_type.lower())
        return namespace


    def _get_uuid(namespace: uuid.UUID, sequence: str) -> uuid.UUID:
        """
        Returns a UUID version 5 for the given namespace and sequence.
        """
        return uuid.uuid5(namespace, sequence)


    def _get_short_uuid(full_uuid: uuid.UUID, sequence: str, length: int = 8) -> str:
        """
        Returns a short, reproducible substring of the UUID for use in accession IDs.
        The substring is selected using a random index seeded by the sequence.
        """
        seed = sum(ord(char) for char in sequence)
        rng = random.Random(seed)
        uuid_str = str(full_uuid).replace("-", "")
        max_start = len(uuid_str) - length
        start = rng.randint(0, max_start)
        short_uuid = uuid_str[start:start + length].upper()
        return short_uuid


    def get_hgrm_accession(sequence: str, accession_type: str) -> tuple[uuid.UUID, str, str]:
        """
        Generates a stable, repeatable accession ID for a given sequence and accession type.
        Args:
            sequence (str): The input sequence for which to generate the UUID.
            accession_type (str): One of 'guide', 'element', 'outer_primer', 'inner_primer'.
        Returns:
            tuple: (uuid.UUID, str, str) - UUID, accession ID, and the original sequence.
        """
        accession_domain = "HGRM"
        accession_type_abbr = _get_accession_type_abbr(accession_type)
        accession_namespace = _get_custom_namespace(accession_type)
        sequence_uuid = _get_uuid(accession_namespace, sequence)
        accession_uuid = _get_short_uuid(sequence_uuid, sequence)
        return sequence_uuid, f"{accession_domain}{accession_type_abbr}{accession_uuid}", sequence
    return (get_hgrm_accession,)


@app.cell
def _(pd):
    guide_metadata_fnp = "metadata/hiPSC-EC/240918_EC_TAPseq_CRISPR.full.tsv"

    guide_metadata_raw = pd.read_csv(guide_metadata_fnp, sep="\t")
    guide_metadata_raw
    return (guide_metadata_raw,)


@app.cell
def _(guide_metadata_raw):
    guide_metadata_raw.shape
    return


@app.cell
def _(mo):
    mo.md(r"""Number of duplicated rows in metadata:""")
    return


@app.cell
def _(guide_metadata_raw):
    guide_metadata_raw.shape[0] - guide_metadata_raw.drop_duplicates().shape[0]
    return


@app.cell
def _(mo):
    mo.md(r"""Number of entries with duplicated pairs (OligoID, CoreOligo) in metadata:""")
    return


@app.cell
def _(guide_metadata_raw):
    guide_metadata_raw.shape[0] - guide_metadata_raw[["OligoID", "CoreOligo"]].drop_duplicates().shape[0]
    return


@app.cell
def _(mo):
    mo.md(r"""Metadata entries with same CoreOligo but different OligoID""")
    return


@app.cell
def _(guide_metadata_raw):
    guide_metadata_raw.groupby('CoreOligo').filter(
        lambda g: g['OligoID'].nunique() == 2
    )
    return


@app.cell
def _(mo):
    mo.md(r"""Deduplicate using CoreOligo and keep the first row is encountered""")
    return


@app.cell
def _(guide_metadata_raw):
    guide_metadata_deduplicated = guide_metadata_raw.drop_duplicates(subset='CoreOligo', keep='first').reset_index(drop=True)

    guide_metadata_deduplicated
    return (guide_metadata_deduplicated,)


@app.cell
def _(guide_metadata_deduplicated):
    guide_bed = guide_metadata_deduplicated[["chr", "start", "end", "name"]]

    req = ["chr", "start", "end", "name"]
    guide_bed = guide_bed.dropna(subset=req).reset_index(drop=True)

    guide_bed = guide_bed.astype({
        "chr":  "string",
        "start":"int64",
        "end":  "int64",
        "name": "string",
    })

    guide_bed
    return (guide_bed,)


@app.cell
def _(guide_bed):
    import pybedtools
    targeted_elements =(
        pybedtools.BedTool.from_dataframe(guide_bed)
        .sort(g="../../Annotations/HG38/GRCh38_EBV.chrom.sizes.no.alt.tsv")
        .merge(d=50, c=4, o='distinct', s=False)
    )
    targeted_elements.to_dataframe()
    return


@app.cell
def _(guide_metadata_deduplicated):
    mask_safe_targeting = guide_metadata_deduplicated["guideSet"] == "safe_targeting"
    mask_negative_control = guide_metadata_deduplicated["guideSet"] == "negative_control"
    mask_targeting = ~((mask_safe_targeting) | (mask_negative_control))
    return mask_negative_control, mask_safe_targeting, mask_targeting


@app.cell
def _(mo):
    mo.md(r"""Number of "safe_targeting" and "negative_control" guides""")
    return


@app.cell
def _(guide_metadata_deduplicated, mask_negative_control, mask_safe_targeting):
    guide_metadata_deduplicated[(mask_safe_targeting) | (mask_negative_control)]["guideSet"].value_counts()
    return


@app.cell
def _(mo):
    mo.md(r"""Number of guides associated with each element""")
    return


@app.cell
def _(guide_metadata_deduplicated):
    # To get the spacer sequence I am going to take the `CoreOligo` and remove the `RightGA` and `LeftGA`
    right_ga = guide_metadata_deduplicated["RightGA"].unique()[0]
    left_ga = guide_metadata_deduplicated["LeftGA"].unique()[0]

    guide_metadata_deduplicated["ordered_sequence"]=guide_metadata_deduplicated["CoreOligo"].str.replace(pat=left_ga, repl="", regex=False).str.replace(pat=right_ga, repl="", regex=False)
    guide_metadata_deduplicated["ordered_sequence"]
    return


@app.cell
def _(guide_metadata_deduplicated):
    guide_metadata_reduced = guide_metadata_deduplicated[["OligoID", "locus", "guideSet","ordered_sequence"]]

    # I am not saving this anymore because incorrect
    #guide_metadata_reduced.to_csv("results/2025-10-01_ec_guide_annotation/metadata_guide_reduced_oligo_stripped.csv", header=True, index=False)
    guide_metadata_reduced
    return (guide_metadata_reduced,)


@app.cell
def _(guide_metadata_reduced):
    guides_to_convert = guide_metadata_reduced[["OligoID", "ordered_sequence"]]
    #
    guides_to_convert
    return (guides_to_convert,)


@app.cell
def _(guides_to_convert):
    # Do we have multiple sequences associated with the same OligoID
    unique_seq_counts = guides_to_convert.groupby('ordered_sequence')['OligoID'].nunique()
    # Get the unique pairs
    guides_to_convert_deduplicated = guides_to_convert[~guides_to_convert.duplicated()]

    # Filter for IDs that have more than one unique sequence
    ids_with_multiple_sequences = unique_seq_counts[unique_seq_counts > 1]

    if not ids_with_multiple_sequences.empty:
        print("\nFound IDs associated with multiple UNIQUE sequences (Potential Data Error):")
        print(ids_with_multiple_sequences)
    else:
        print("\nAll IDs map to a single unique sequence.")
    return


@app.cell
def _(mo):
    mo.md(
        r"""
    ### Inconsistencies in guide file
    The guide designed file contains a lot of duplicated entries.
    The same sequence oligo is associated with multiple IDs.
    Completely disregard the file and use the `feature_reference_file.csv` that is passed to `cellranger`.
    """
    )
    return


@app.cell
def _(mo):
    mo.md(r"""# Generating FASTQ files for alignment""")
    return


@app.cell
def _(pd):
    # I am going to use the `guide_metadata_reduced` dataframe to plot annotations
    feature_reference = pd.read_csv("metadata/hiPSC-EC/iPSC-EC_feature-reference.csv")
    feature_reference
    return (feature_reference,)


@app.cell
def _(guide_metadata_reduced):
    guide_metadata_reduced
    return


@app.cell
def _(feature_reference, get_hgrm_accession, guide_metadata_reduced):
    def _():
        ### CHECK AND REMOVE
        _tmp_for_merging = guide_metadata_reduced.rename(columns={"OligoID":"name", "ordered_sequence":"sequence" })

        guide_metadata_deduplicated = (
            feature_reference[["name","sequence"]]
            .merge(_tmp_for_merging, how="left", on="name")
            .rename(columns={"sequence_x":"spacer_with_G", "sequence_y":"spacer" })
        )



        # Add HGRM ID by applying the function get_hgrm_accession(sequence1, "guide") to all the sequences in the spacer column

        guide_metadata_deduplicated[['uuid', 'accession_id']] = guide_metadata_deduplicated.apply(
            lambda row: get_hgrm_accession(row[4], "guide")[:2],
            axis=1,
            result_type='expand'
        )

        guide_metadata_deduplicated.to_csv("results/2025-10-01_ec_guide_annotation/metadata_guides_with_accession.csv", index=False)
        return guide_metadata_deduplicated


    _()
    return


@app.cell
def _(alt, guide_metadata_deduplicated, mask_targeting, pd):
    _mask_safe_targeting = guide_metadata_deduplicated["guideSet"] == "safe_targeting"
    _mask_negative_control = guide_metadata_deduplicated["guideSet"] == "negative_control"
    _mask_targeting = ~((_mask_safe_targeting) | (_mask_negative_control))

    guide_metadata_deduplicated[(_mask_safe_targeting) | (_mask_negative_control)]["guideSet"].value_counts()

    guide_set_counts = guide_metadata_deduplicated[mask_targeting]["guideSet"].value_counts()

    # 1. Data Preparation: Convert the Series into a DataFrame suitable for Altair
    plot_df = guide_set_counts.reset_index(name='Count')
    plot_df.columns = ['GuideSet_ID', 'Count'] # Renaming columns

    # 2. Calculate Mean and Median for annotation rules
    mean_val = plot_df['Count'].mean()
    median_val = plot_df['Count'].median()

    # Create a DataFrame for the annotation rules (mean and median lines)
    stats_df = pd.DataFrame({
        'value': [mean_val, median_val],
        'label': [f'Mean: {mean_val:.0f}', f'Median: {median_val:.0f}']
    })


    histogram = alt.Chart(plot_df).encode(
        x=alt.X('Count', 
                bin=alt.Bin(maxbins=50), 
                title="Guide Set Size (Count)", 
                # --- FIX: Use tickBand="center" to align ticks ---
                axis=alt.Axis(
                    tickBand="center" 
                )
            ),
        y=alt.Y('count()', title="Frequency of Sizes")
    ).mark_bar()

    # 2. Annotation Rules (Mean and Median Lines - remains unchanged)
    annotations = alt.Chart(stats_df).mark_rule(size=2).encode(
        x='value',
        color=alt.Color('label', scale=alt.Scale(range=['#1B9E77', '#D95F02'])),
        tooltip=[alt.Tooltip('label'), alt.Tooltip('value', format=',.0f', title='Value')]
    )

    # 3. Layer and Save the Chart
    chart = (histogram + annotations).properties(
        title="Distribution of Guide Set Sizes"
    )

    chart.save('plots/distribution_designed_ngrnas_per_element.png', scale_factor=4)

    chart
    return


@app.cell
def _(feature_reference):
    feature_reference.columns
    return


@app.cell
def _(feature_reference):
    # Write the FASTQ file using I for scores
    fastq_records_series = feature_reference[["id", "sequence"]].apply(
            lambda row: (
                # Line 1: Identifier (@Name)
                f"@{row.iloc[0]}\n" 
                # Line 2: Sequence
                f"{row.iloc[1][1:]}\n"    
                # Line 3: Separator
                f"+\n"                
                # Line 4: Quality Scores ('I' * length)
                f"{'I' * len(row.iloc[1][1:])}" 
            ),
            axis=1
        )
    # Create the string to write
    final_content = '\n'.join(fastq_records_series.tolist())

    # Write the final string to the output file
    with open("results/2025-10-01_ec_guide_annotation/ec_differentiation.fastq", 'w') as f:
        f.write(final_content + '\n')
    return


@app.cell
def _(feature_reference):
    def _():
        # Write the FASTQ file using I for scores. Adding the PAM
        fastq_records_series = feature_reference[["id", "sequence"]].apply(
                lambda row: (
                    # Line 1: Identifier (@Name)
                    f"@{row.iloc[0]}\n" 
                    # Line 2: Sequence
                    f"{row.iloc[1][1:]}NGG\n"    
                    # Line 3: Separator
                    f"+\n"                
                    # Line 4: Quality Scores ('I' * length)
                    f"{'I' * (int(len(row.iloc[1][1:]))+3)}" 
                ),
                axis=1
            )
        # Create the string to write
        final_content = '\n'.join(fastq_records_series.tolist())

        # Write the final string to the output file
        with open("results/2025-10-01_ec_guide_annotation/ec_differentiation_with_pam.fastq", 'w') as f:
            return f.write(final_content + '\n')


    _()
    return


@app.cell
def _(feature_reference):
    def _():
        # Write the FASTQ file using I for scores. Adding the PAM
        fastq_records_series = feature_reference[["id", "sequence"]].apply(
                lambda row: (
                    # Line 1: Identifier (@Name)
                    f"@{row.iloc[0]}\n" 
                    # Line 2: Sequence
                    f"{row.iloc[1]}NGG\n"    
                    # Line 3: Separator
                    f"+\n"                
                    # Line 4: Quality Scores ('I' * length)
                    f"{'I' * (int(len(row.iloc[1]))+3)}" 
                ),
                axis=1
            )
        # Create the string to write
        final_content = '\n'.join(fastq_records_series.tolist())

        # Write the final string to the output file
        with open("results/2025-10-01_ec_guide_annotation/ec_differentiation_with_pam_and_leading_g.fastq", 'w') as f:
            return f.write(final_content + '\n')


    _()
    return


if __name__ == "__main__":
    app.run()
