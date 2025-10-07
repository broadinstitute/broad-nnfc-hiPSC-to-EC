import marimo

__generated_with = "0.15.2"
app = marimo.App()


@app.cell
def _():
    import marimo as mo
    return


@app.cell
def _():
    import altair as alt
    import pandas as pd
    import numpy as np
    import csv
    import pybedtools
    import pysam
    return alt, csv, np, pd, pybedtools


@app.cell
def _():
    guide_merge_threshold = 100 # Merge together guides closer than dist_threshold
    promoter_distance = 500 # Define the promoter as +/- 500bp from TSS
    regulatory_region_distance = 1_000_000 # Extend +/- 1Mb around the TSS
    return guide_merge_threshold, promoter_distance, regulatory_region_distance


@app.cell
def _():
    # Guides not aligning to hg38
    unmapped_guide_file = "results/2025-10-01_ec_guide_annotation/ec_differentiation_with_pam_and_leading_g_hg38_unmapped.tsv"
    # Guides that are mapped after aligning to hg38
    mapped_guide_file = "results/2025-10-01_ec_guide_annotation/ec_differentiation_with_pam_and_leading_g_hg38_filtered.bed"
    # Chromosome size file
    genome_file = "../../Annotations/HG38/GRCh38_EBV.chrom.sizes.no.alt.tsv"
    # TSS annotation file from GENCODE v44
    tss_file = "results/2025-09-15/intermediates/gencode_tss.bed"
    # Guide metadata file
    guide_metadata_file = "metadata/hiPSC-EC/metadata_guides_with_accession.csv"
    # Gene metadata file
    gene_metadata_file = "metadata/hiPSC-EC/metadata_genes_corrected_names.csv"
    return (
        gene_metadata_file,
        genome_file,
        guide_metadata_file,
        mapped_guide_file,
        tss_file,
        unmapped_guide_file,
    )


@app.cell
def _(gene_metadata_file, pd):
    # Get the list of dial-out genes
    readout_genes = pd.read_csv(
        gene_metadata_file
    )
    readout_genes =(
        readout_genes["intended_target_gene_symbol"]
        .drop_duplicates()
        .reset_index(drop=True)
    )
    readout_genes = readout_genes.to_frame()
    readout_genes
    return (readout_genes,)


@app.cell
def _(genome_file, pybedtools, tss_file):
    # Loading the gene TSS annotation
    tss_bed = pybedtools.BedTool(tss_file).sort(g=genome_file)
    tss_bed.head()
    return (tss_bed,)


@app.cell
def _(mapped_guide_file, pybedtools):
    pybedtools.BedTool(mapped_guide_file).tail()
    return


@app.cell
def _(guide_metadata_file, pd):
    # Read the guide metadata
    guide_metadata = pd.read_csv(guide_metadata_file).rename(columns={"name":"guide_id"})
    # Mapping guide_id and accession_id
    guide_accession_map = guide_metadata[['guide_id', 'accession_id']]

    guide_metadata.head()
    return guide_accession_map, guide_metadata


@app.cell
def _(guide_accession_map, pd, unmapped_guide_file):
    # Read the unmapped guides file
    unmapped_guides = pd.read_csv(unmapped_guide_file, sep='\t')

    # Perform the left merge to add the 'accession_id' column
    unmapped_guides = unmapped_guides.merge(
        guide_accession_map,
        on='guide_id',
        how='left'
    )
    unmapped_guides.head()
    return


@app.cell
def _(mapped_guide_file, pybedtools):
    pybedtools.BedTool(mapped_guide_file).tail()
    return


@app.cell
def _(genome_file, guide_merge_threshold, mapped_guide_file, pybedtools):
    targeted_elements =(
        pybedtools.BedTool(mapped_guide_file)
        .sort(g=genome_file)
        .merge(d=guide_merge_threshold, c=4, o='distinct', s=False)
    )
    len(targeted_elements)
    return (targeted_elements,)


@app.cell
def _(csv, guide_accession_map, targeted_elements):
    grna_to_element = targeted_elements.to_dataframe()
    # Save the gRNA to element mapping
    grna_to_element['grna_target'] = (
        grna_to_element['chrom']
            .str.cat(
                (grna_to_element['start'] + 1).astype('int64').astype(str), sep=':'
            ).str.cat(
                (grna_to_element['end']).astype('int64').astype(str), sep='-'
            )
    )
    exploded_grna_to_element = grna_to_element.assign(
        name=grna_to_element['name'].str.split(',')
    ).explode('name')

    final_mapping = exploded_grna_to_element[['name', 'grna_target']].drop_duplicates().rename(columns={'name': 'grna_id'})

    # Use accession_id instead of grna_id
    final_mapping = final_mapping.merge(
        guide_accession_map,
        left_on='grna_id',
        right_on='guide_id',
        how='left'
    ).drop(columns=['grna_id', 'guide_id']).rename(columns={'accession_id': 'grna_id'})[['grna_id', 'grna_target']]
    # write to a TSV file
    final_mapping.to_csv(
        "results/2025-10-02_ec_create_sceptre_input/grna_to_element_mapping_hg38.tsv",
        sep='\t',
        index=False,
        header=True,
        quoting=csv.QUOTE_NONE
    )
    final_mapping
    return (final_mapping,)


@app.cell
def _(alt, final_mapping, pd):
    final_mapping_counts = final_mapping["grna_target"].value_counts()

    # 1. Data Preparation: Convert the Series into a DataFrame suitable for Altair
    plot_df = final_mapping_counts.reset_index(name='Count')
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

    chart.save('plots/distribution_designed_ngrnas_per_element_after_alignment.png', scale_factor=4)

    chart
    return


@app.cell
def _(guide_metadata):
    # Count the total occurrences of 'locus'
    total_guides_annotation = (
        guide_metadata['locus']
            .value_counts()
            .to_frame('total_guide_count')
    )

    total_guides_annotation
    return (total_guides_annotation,)


@app.cell
def _(targeted_elements):
    targeted_elements.to_dataframe()
    return


@app.cell
def _(targeted_elements):
    element_guide_counts = (
        targeted_elements
        .to_dataframe(usecols=[3])
    )
    element_guide_counts['number_of_guides'] = element_guide_counts['name'].str.split(',').apply(len)
    element_guide_counts
    return


@app.cell
def _(alt, targeted_elements):
    def _():
        # Distribution of number of guides per element
        element_guide_counts = (
            targeted_elements
            .to_dataframe(usecols=[3])
        )
        element_guide_counts['number_of_guides'] = element_guide_counts['name'].str.split(',').apply(len)

        num_elements = len(element_guide_counts)

        chart = (
            alt.Chart(element_guide_counts)
            .mark_bar()
            .encode(
                x=alt.X('number_of_guides', bin=alt.Bin(maxbins=30), title='Number of Guides per Element'),
                y=alt.Y('count()', title='Element Count')
            )
            .properties(
                title=alt.Title(text="Distribution of Number of Guides per Targeted Element",subtitle=f"(n={num_elements})")
            )
        )
        return chart


    _()
    return


@app.cell
def _(tss_bed):
    gene_id_name_map =(
        tss_bed
            .to_dataframe(
                header=None,
                usecols=[6,8]
            )
            .drop_duplicates()
            .reset_index(drop=True)
            .rename(
                columns={"thickStart": 'gene_id', "itemRgb": 'gene_symbol'},
            )
    )
    gene_id_name_map.head()
    return (gene_id_name_map,)


@app.cell
def _(gene_id_name_map, readout_genes):
    m = (
        gene_id_name_map
         .drop_duplicates('gene_symbol')          # handle duplicates
         .set_index('gene_symbol')['gene_id']
    )


    # optional: normalize case/whitespace
    readout_genes['intended_target_gene_symbol'] = readout_genes['intended_target_gene_symbol'].str.strip()
    m.index = m.index.str.strip()

    # vectorized map; missing symbols -> NaN (no KeyError)
    readout_genes['gene_id'] = (
        readout_genes['intended_target_gene_symbol']
            .map(m)
    )
    readout_genes.head(), len(readout_genes)
    return


@app.cell
def _(genome_file, regulatory_region_distance, tss_bed):
    # Create regulatory regions as +/- 1Mb around TSS
    # Use regulatory_region_distance at the top of the notebook to set this value
    regulatory_region_pybed = tss_bed.slop(
        g=genome_file,
        s=True,
        l=regulatory_region_distance,
        r=regulatory_region_distance
    ).sort(g=genome_file)
    return (regulatory_region_pybed,)


@app.cell
def _(genome_file, promoter_distance, tss_bed):
    # Create promoter regions as +/- 500bp around TSS
    # Use promoter_distance at the top of the notebook to set this value
    promoter_region_pybed = tss_bed.slop(
        g=genome_file,
        s=True,
        l=promoter_distance,
        r=promoter_distance
    ).sort(g=genome_file)
    return (promoter_region_pybed,)


@app.cell
def _(genome_file, regulatory_region_pybed, targeted_elements):
    # Find all elements in the regulatory region for every gene (not just readout genes)
    # Require 50% of the element to overlap
    e_g_pairs = (
        regulatory_region_pybed
            .intersect(
                targeted_elements,
                g=genome_file,
                wa=True,
                wb=True,
                F=0.5,
            )
    )
    e_g_pairs.head()
    return (e_g_pairs,)


@app.cell
def _(e_g_pairs):
    # Format the results in a dataframe
    discovery_pairs = (
        e_g_pairs
            .to_dataframe(header=None, usecols=[6,9,10,11])
            .rename(
                columns={6: 'gene_id',
                         9: 'element_chrom',
                         10: 'element_start',
                         11: 'element_end'}
            )
    )
    # Create the element name as chrom:start-end
    discovery_pairs['grna_target'] = (
        discovery_pairs['element_chrom']
          .str.cat(
              (discovery_pairs['element_start'] + 1).astype('int64').astype(str), sep=':'
          )
          .str.cat(
              (discovery_pairs['element_end']).astype('int64').astype(str), sep='-'
          )
    )
    # 
    discovery_pairs = (
        discovery_pairs
            .drop(columns=['element_chrom','element_start','element_end'])
            .drop_duplicates()
            .reset_index(drop=True)
    )
    discovery_pairs.head()
    return (discovery_pairs,)


@app.cell
def _(e_g_pairs, genome_file, promoter_region_pybed, targeted_elements):
    # Find all elements in the regulatory region for every gene (not just readout genes)
    # Require 50% of the element to overlap
    e_tss_pairs = (
        promoter_region_pybed
            .intersect(
                targeted_elements,
                g=genome_file,
                wa=True,
                wb=True,
                F=0.5,
            )
    )
    e_g_pairs.head()
    return (e_tss_pairs,)


@app.cell
def _(e_tss_pairs):
    # Format the results in a dataframe
    positive_pairs = (
        e_tss_pairs
            .to_dataframe(header=None, usecols=[6,9,10,11])
            .rename(
                columns={6: 'gene_id',
                         9: 'element_chrom',
                         10: 'element_start',
                         11: 'element_end'}
            )
    )
    # Create the element name as chrom:start-end
    positive_pairs['grna_target'] = (
        positive_pairs['element_chrom']
          .str.cat(
              (positive_pairs['element_start'] + 1).astype('int64').astype(str), sep=':'
          )
          .str.cat(
              (positive_pairs['element_end']).astype('int64').astype(str), sep='-'
          )
    )
    # 
    positive_pairs = (
        positive_pairs
            .drop(columns=['element_chrom','element_start','element_end'])
            .drop_duplicates()
            .reset_index(drop=True)
    )
    positive_pairs.head()
    return (positive_pairs,)


@app.cell
def _(discovery_pairs, positive_pairs):
    # Remove all the pairs that overlaps with positive pairs
    keys = ['grna_target', 'gene_id']

    discovery_pairs_no_pos = (
        discovery_pairs
            .merge(
                positive_pairs[keys].drop_duplicates(),
                on=keys,
                how='left',
                indicator=True,
            )
            .query('_merge == "left_only"')
            .drop(columns='_merge')
            .reset_index(drop=True)
    )
    discovery_pairs_no_pos.head(), len(discovery_pairs_no_pos)
    return discovery_pairs_no_pos, keys


@app.cell
def _(csv, discovery_pairs_no_pos, keys, positive_pairs, readout_genes):
    # Keep only readout genes
    discovery_pairs_final = (
        discovery_pairs_no_pos
            .merge(
                readout_genes[['gene_id']],
                on='gene_id'
            )[keys]
            .rename(columns={"gene_id": "response_id"})
    )

    positive_pairs_final = (
        positive_pairs
            .merge(
                readout_genes[['gene_id']],
                on='gene_id'
            )[keys]
            .rename(columns={"gene_id": "response_id"})
    )

    discovery_pairs_final.to_csv(
        "results/2025-09-15/discovery_pairs_no_pos_hg38_readout_genes.tsv",
        sep='\t',
        index=False,
        header=False,
        quoting=csv.QUOTE_NONE
    )

    positive_pairs_final.to_csv(
        "results/2025-09-15/positive_pairs_hg38_readout_genes.tsv",
        sep='\t',
        index=False,
        header=False,
        quoting=csv.QUOTE_NONE
    )
    return discovery_pairs_final, positive_pairs_final


@app.cell
def _(positive_pairs_final):
    positive_pairs_final
    return


@app.cell
def _(discovery_pairs_final, positive_pairs_final, targeted_elements):
    # Get all grna_target values from both filtered dataframes
    used_grna_targets = (
        set(
            discovery_pairs_final["grna_target"]
        ).union(
            positive_pairs_final["grna_target"]
        )
    )

    targeted_elements_df = targeted_elements.to_dataframe()

    targeted_elements_df['grna_target'] = (
        targeted_elements_df['chrom']
            .str.cat(
                (targeted_elements_df['start'] + 1).astype('int64').astype(str), sep=':'
            )
            .str.cat(
                (targeted_elements_df['end']).astype('int64').astype(str), sep='-'
            )
    )
    # Find grna_targets in targeted_elements_df not present in either dataframe
    used_grna_mask = targeted_elements_df['grna_target'].isin(used_grna_targets)

    unused_elements_df = targeted_elements_df[~used_grna_mask]

    all_unused_guides = unused_elements_df['name'].str.split(',').explode().unique()
    all_unused_guides
    return (all_unused_guides,)


@app.cell
def _(all_unused_guides):
    len(all_unused_guides)
    return


@app.cell
def _(all_unused_guides, guide_metadata, total_guides_annotation):
    len(all_unused_guides)
    # Count the occurrences of 'locus' for guides with bad reads
    unused_guides_mask = guide_metadata["guide_id"].isin(list(all_unused_guides))

    unused_guides_annotation = (
        guide_metadata[unused_guides_mask]['locus']
            .value_counts()
            .to_frame('safe_targeting_count')
    )

    guide_count_per_element = total_guides_annotation.join(unused_guides_annotation, how='left').fillna(0).astype(int)
    guide_count_per_element
    return


@app.cell
def _(genome_file, targeted_elements, tss_bed):
    closest_genes = targeted_elements.closest(tss_bed, d=True, g=genome_file)
    closest_genes.head()
    return (closest_genes,)


@app.cell
def _(closest_genes, np):
    cols = ['element_chrom','element_start','element_end', 'element_name', 'tss_chrom','tss_start','tss_end','transcript_id','score','transcript_strand', 'gene_id', 'gene_type', 'gene_symbol', 'dist'
    ]
    df = closest_genes.to_dataframe(names=cols)

    # Use region midpoint vs TSS position (TSS is 1 bp; start == position in 0-based BED)
    df['element_mid'] = (df['element_start'] + df['element_end']) // 2
    df['tss_pos'] = df['tss_start']

    # signed distance relative to transcription direction
    sign = df['transcript_strand'].map({'+': 1, '-': -1})
    df['signed_dist'] = (df['element_mid'] - df['tss_pos']) * sign

    # classify side (overlaps -> special case)
    df['tss_side'] = np.where(
        df['dist'].eq(0), 'overlap',
        np.where(df['signed_dist'] < 0, 'upstream', 'downstream')
    )

    # optional: keep only what you need
    out = df[['element_chrom','element_start','element_end','tss_start','transcript_strand','gene_id','dist','signed_dist','tss_side']]
    out.tss_side.value_counts()
    return


if __name__ == "__main__":
    app.run()
