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
    return alt, csv, np, pd, pybedtools, pysam


@app.cell(hide_code=True)
def _(pysam):
    def process_alignments_and_track_unmapped(sam_file_path, bed_output_path):
        """
        Processes a SAM file in a single pass to create a BED file of good alignments
        and return a set of names for unmapped/non-canonical reads.

        Args:
            sam_file_path (str): Path to the input SAM/BAM file.
            bed_output_path (str): Path for the output BED file.

        Returns:
            set: A set of unique read names considered unmapped for analysis.
        """
        # Define canonical chromosomes
        canonical_chroms = {f'chr{i}' for i in range(1, 23)}
        canonical_chroms.update(['chrX', 'chrY', 'chrM'])

        read_info = {}
        reads_written = 0

        try:
            with pysam.AlignmentFile(sam_file_path, "r") as samfile, \
                 open(bed_output_path, "w") as bedfile:

                for read in samfile:
                    read_name = read.query_name

                    # Store info for later tracking of unmapped/non-canonical reads
                    if read_name not in read_info:
                        read_info[read_name] = {'mapped_chroms': set(), 'is_unmapped_flag': False}

                    # Task 1: Write "good" alignments to the BED file
                    is_primary = not read.is_unmapped and not read.is_secondary and not read.is_supplementary
                    is_canonical = read.reference_name in canonical_chroms if not read.is_unmapped else False

                    if is_primary and is_canonical:
                        # Extract and write to BED file
                        bed_line = f"{read.reference_name}\t{read.reference_start}\t{read.reference_end}\t" \
                                   f"{read_name}\t{read.mapping_quality}\t{'+' if not read.is_reverse else '-'}\n"
                        bedfile.write(bed_line)
                        reads_written += 1

                    # Task 2: Keep track of all alignments for a read
                    if read.is_unmapped:
                        read_info[read_name]['is_unmapped_flag'] = True
                    else:
                        read_info[read_name]['mapped_chroms'].add(read.reference_name)

        except FileNotFoundError:
            print(f"Error: The file {sam_file_path} was not found.")
            return set()

        # Final pass to identify "bad" reads based on collected info
        unmapped_read_names = set()
        for read_name, info in read_info.items():
            # A read is "bad" if it's unmapped OR all its alignments are non-canonical
            if info['is_unmapped_flag'] or not info['mapped_chroms'].intersection(canonical_chroms):
                unmapped_read_names.add(read_name)

        print(f"Wrote {reads_written} canonical primary alignments to {bed_output_path}.")
        print(f"Identified {len(unmapped_read_names)} reads as unmapped or non-canonical.")

        return unmapped_read_names
    return


@app.cell
def _():
    guide_merge_threshold = 50 # Merge together guides closer than dist_threshold
    promoter_distance = 500 # Define the promoter as +/- 500bp from TSS
    regulatory_region_distance = 1_000_000 # Extend +/- 1Mb around the TSS
    return guide_merge_threshold, promoter_distance, regulatory_region_distance


@app.cell
def _():
    bed_file = "metadata/hiPSC-HPC/guides_original_positions_hg19_lifted_to_hg38.bed"
    genome_file = "../../Annotations/HG38/GRCh38_EBV.chrom.sizes.no.alt.tsv"
    tss_file = "results/2025-09-15/gencode_tss.bed"
    guide_metadata_file = "metadata/hiPSC-HPC/All_guides_with_addnl_controls_altTSS.csv"
    gene_metadata_file = "metadata/hiPSC-HPC/metadata_genes.csv"
    return (
        bed_file,
        gene_metadata_file,
        genome_file,
        guide_metadata_file,
        tss_file,
    )


@app.cell
def _(
    gene_metadata_file,
    genome_file,
    guide_metadata_file,
    pd,
    pybedtools,
    tss_file,
):
    readout_genes = pd.read_csv(
        gene_metadata_file,
        header=1,
        usecols=[6],
        names=['gene_symbol'],
    )

    tss_bed = pybedtools.BedTool(tss_file).sort(g=genome_file)

    guide_metadata = pd.read_csv(guide_metadata_file).set_index("OligoID")

    guide_metadata.head()
    return guide_metadata, readout_genes, tss_bed


@app.cell
def _(bed_file, pybedtools):
    pybedtools.BedTool(bed_file).tail()
    return


@app.cell
def _(bed_file, genome_file, guide_merge_threshold, pybedtools):
    targeted_elements =(
        pybedtools.BedTool(bed_file)
        .sort(g=genome_file)
        .merge(d=guide_merge_threshold, c=4, o='distinct')
    )
    len(targeted_elements)
    return (targeted_elements,)


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
def _(alt, targeted_elements):
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
            x=alt.X('number_of_guides', bin=alt.Bin(maxbins=50), title='Number of Guides per Element'),
            y=alt.Y('count()', title='Element Count')
        )
        .properties(
            title=alt.Title(text="Distribution of Number of Guides per Targeted Element",subtitle=f"(n={num_elements})")
        )
    )
    chart
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
    readout_genes['gene_symbol'] = readout_genes['gene_symbol'].str.strip()
    m.index = m.index.str.strip()

    # vectorized map; missing symbols -> NaN (no KeyError)
    readout_genes['gene_id'] = (
        readout_genes['gene_symbol']
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
def _(all_unused_guides, guide_metadata, total_guides_annotation):
    len(all_unused_guides)
    # Count the occurrences of 'locus' for guides with bad reads
    unused_guides_mask = guide_metadata.index.isin(all_unused_guides)

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
