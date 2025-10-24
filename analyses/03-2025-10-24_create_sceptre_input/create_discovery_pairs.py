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
    import numpy as np
    import csv
    import pybedtools
    import pysam
    return alt, csv, np, pd, pybedtools


@app.cell
def _():
    from scripts.python.create_accession_id import get_accession
    return (get_accession,)


@app.cell
def _():
    guide_merge_threshold = 50 # Merge together guides closer than dist_threshold
    promoter_distance = 500 # Define the promoter as +/- 500bp from TSS
    regulatory_region_distance = 1_000_000 # Extend +/- 1Mb around the TSS
    return guide_merge_threshold, promoter_distance, regulatory_region_distance


@app.cell
def _():
    # Guides not aligning to hg38
    unmapped_guide_file = "results/02-2025-10-21_align_guides/unmapped.tsv"
    # Not valid alignments
    discarded_alignments_file = "results/02-2025-10-21_align_guides/discarded.tsv"
    # Guides that are mapped after aligning to hg38
    mapped_guide_file = "results/02-2025-10-21_align_guides/ec_differentiation_with_pam_and_leading_g_hg38_filtered.bed"
    # Chromosome size file
    genome_file = "../../Annotations/HG38/GRCh38_EBV.chrom.sizes.no.alt.tsv"
    # TSS annotation file from GENCODE v44
    tss_file = "results/2025-09-15/intermediates/gencode_tss.bed"
    # Guide metadata file
    cellranger_feature_reference_file = "metadata/iPSC-EC_feature-reference.csv"
    # DEPRECATED: Raw metadata guides. File unusable.
    #legacy_metadata_file = "metadata/240918_EC_TAPseq_CRISPR.full.tsv"
    # Gene metadata file
    gene_metadata_file = "metadata/metadata_genes_final.csv"
    return (
        cellranger_feature_reference_file,
        discarded_alignments_file,
        gene_metadata_file,
        genome_file,
        mapped_guide_file,
        tss_file,
        unmapped_guide_file,
    )


@app.cell
def _(mo):
    mo.md(
        r"""
    ## Load guide sequences metadata
    I am going to use the feature reference file passed as input to cellranger to map the old guide IDs to sequences. The legacy metadata has inconsistent accessions and contains redundancy and multiple accessions for the same sequence.
    I am assigning new accessions that are spacer based. Same spacer will give the same accession.
    """
    )
    return


@app.cell
def _(cellranger_feature_reference_file, csv, get_accession, mo, pd):
    cellranger_feature_reference = pd.read_csv(cellranger_feature_reference_file, sep=',')
    cellranger_feature_reference['sequence'] = cellranger_feature_reference['sequence'].str.upper().str[1:]

    message = (
        f"Loaded {len(cellranger_feature_reference)} guide sequences "
        "from CellRanger feature reference."
    )
    mo.md(message)
    guide_metadata_df = cellranger_feature_reference[['sequence', 'id']].rename(columns={'sequence':'spacer', 'id': 'aliases' })

    # The function get_accession_id returns three values: sequence_uuid, guide_id, and the guide sequence.
    ACCESSION_TYPE = 'GUIDE'
    tmp = guide_metadata_df['spacer'].map(
        lambda s: get_accession(s, accession_type=ACCESSION_TYPE)
    ).apply(pd.Series)

    tmp.columns = ['guide_uuid', 'guide_id', 'guide_sequence']  # adjust if your function returns other order

    # Attach just the UUID/ID (drop the returned sequence unless you want it)
    guide_metadata_df = pd.concat(
        [tmp[['guide_id', 'guide_uuid']], guide_metadata_df], axis=1
    )

    # Save the final guide metadata with consistent accessions
    guide_metadata_df.to_csv(
        "results/03-2025-10-24_create_sceptre_input/guide_metadata_with_new_accessions.tsv",
        sep='\t',
        index=False,
        header=True,
        quoting=csv.QUOTE_NONE
    )

    guide_metadata_df.head()
    return (guide_metadata_df,)


@app.cell
def _(guide_metadata_df, pd):
    # Create a map from aliases to guide_id
    alias_to_guide_id = dict(zip(guide_metadata_df['aliases'], guide_metadata_df['guide_id']))
    pd.Series(alias_to_guide_id).head(10)
    return (alias_to_guide_id,)


@app.cell
def _(mo):
    mo.md(
        r"""
    ## Load genes metadata
    Loading the list of readout genes and the TSS for all the genes in GENCODE v44.
    """
    )
    return


@app.cell
def _(gene_metadata_file, pd):
    # Save # Get the list of readout genes
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
def _(mo):
    mo.md(
        r"""
    ## Load the alignment results
    The alignment and filtering steps create three files:

     - `ec_differentiation_with_pam_and_leading_g_hg38_filtered.bed`
     - `discarded.tsv`
     - `unmapped.tsv`

    The `discarded.tsv` file contains all the alignments that are not valid. All the guides that don't have a valid alignment need to be added to the unmapped guides list.
    """
    )
    return


@app.cell
def _(mo):
    mo.md(r"""Load the BED file with the results of the alignment. Convert the names to accessions and create the pybedtools object.""")
    return


@app.cell
def _(alias_to_guide_id, mapped_guide_file, pd, pybedtools):
    # Read the mapped guides BED file
    mapped_guides_df = pd.read_csv(mapped_guide_file, sep='\t', header=None)
    mapped_guides_df.columns = ['chrom', 'start', 'end', 'name', 'score', 'strand', "NM", "AS"]
    # Keep primary chromosomes only
    primary_chromosomes = [f'chr{i}' for i in range(1,23)] + ['chrX', 'chrY', 'chrM']
    mapped_guides_df = mapped_guides_df[mapped_guides_df['chrom'].isin(primary_chromosomes)]
    # Convert names to accessions
    # Because it is a BED file I am not changing the header to guide_id
    mapped_guides_df["name"] = mapped_guides_df["name"].map(alias_to_guide_id)
    mapped_guides_df.head()

    mapped_guides = pybedtools.BedTool.from_dataframe(mapped_guides_df)
    mapped_guides.head()
    return mapped_guides, mapped_guides_df


@app.cell
def _(mo):
    mo.md(r"""Load all the discarded alignments and check which guides are not present at least once in the mapped file.""")
    return


@app.cell
def _(alias_to_guide_id, discarded_alignments_file, pd):
    # Load the discarded alignments
    discarded_alignments = pd.read_csv(discarded_alignments_file, sep='\t')
    discarded_alignments.head()
    # Map to guide IDs
    discarded_alignments['guide_id'] = discarded_alignments['read_name'].map(alias_to_guide_id)
    discarded_alignments.head()
    return (discarded_alignments,)


@app.cell
def _(mo):
    mo.md(r"""Find all the guides that don't have a valid alignment.""")
    return


@app.cell
def _(discarded_alignments, mapped_guides_df):
    kept_guides = set(mapped_guides_df['name'])
    # Guides that are mapped but without a single valid alignment are considered unmapped
    discarded_guides = (set(discarded_alignments['guide_id']) - kept_guides)
    len(discarded_guides)
    return discarded_guides, kept_guides


@app.cell
def _(alias_to_guide_id, pd, unmapped_guide_file):
    # Read the unmapped file
    unmapped_guides = pd.read_csv(unmapped_guide_file, sep='\t')
    unmapped_guides
    # Map these to accession IDs
    unmapped_guides['guide_id'] = unmapped_guides['read_name'].map(alias_to_guide_id)

    # Add a reason column
    unmapped_guides['reason'] = 'unmapped'
    unmapped_guides = unmapped_guides[['guide_id', 'reason']]
    unmapped_guides.head()
    return (unmapped_guides,)


@app.cell
def _(discarded_guides, kept_guides, pd, unmapped_guides):
    all_guides = pd.concat([
        pd.DataFrame(list(kept_guides), columns=['guide_id']).assign(reason='kept'),
        pd.DataFrame(list(discarded_guides), columns=['guide_id']).assign(reason='discarded'),
        unmapped_guides
    ])
    all_guides.groupby('reason')['guide_id'].count()
    return (all_guides,)


@app.cell
def _(all_guides):
    all_guides.shape
    return


@app.cell
def _(mo):
    mo.md(
        r"""
    ## Create list of assayed elements
    Now that we have the coordinates of the guides, we can merge the ones close to each other to create the list of elements being targeted by the same guides.
    Parameters used in the analysis:

     - guide_merge_threshold = 50 # Merge together guides closer than dist_threshold
     - promoter_distance = 500 # Define the promoter as +/- 500bp from TSS
     - regulatory_region_distance = 1_000_000 # Extend +/- 1Mb around the TSS

    """
    )
    return


@app.cell
def _(genome_file, guide_merge_threshold, mapped_guides, mo):
    targeted_elements =(
        mapped_guides
        .sort(g=genome_file)
        .merge(d=guide_merge_threshold, c=4, o='distinct', s=False)
    )

    mo.md(text=f"Number of elements targeted: {len(targeted_elements)}")
    return (targeted_elements,)


@app.cell
def _(csv, targeted_elements):
    targeted_elements_df = targeted_elements.to_dataframe()
    # Save the gRNA to element mapping
    targeted_elements_df['grna_target'] = (
        targeted_elements_df['chrom']
            .str.cat(
                (targeted_elements_df['start'] + 1).astype('int64').astype(str), sep=':'
            ).str.cat(
                (targeted_elements_df['end']).astype('int64').astype(str), sep='-'
            )
    )
    targeted_elements_df['size'] = targeted_elements_df['end'] - targeted_elements_df['start']

    targeted_elements_df['number_of_grna'] = targeted_elements_df['name'].str.split(",").str.len()

    targeted_elements_df = targeted_elements_df[["chrom","start", "end", "grna_target", "size", "name","number_of_grna"]].rename(columns = {"name":"grna_id_list"})

    targeted_elements_df.to_csv("results/03-2025-10-24_create_sceptre_input/targeted_elements.tsv",
            sep="\t",
            index=False,
            header=True,
            quoting=csv.QUOTE_NONE,
           )
    targeted_elements_df.head()
    return (targeted_elements_df,)


@app.cell
def _(alt, targeted_elements_df):
    # Plot distribution of grna per targeted element
    alt.Chart(targeted_elements_df).encode(
        x=alt.X('number_of_grna', 
                bin=alt.Bin(maxbins=50), 
                title="Number of guide RNAs per targeted element (Count)", 
                axis=alt.Axis(
                    tickBand="center" 
                )
            ),
        y=alt.Y('count()', title="Count")
    ).mark_bar()
    return


@app.cell
def _(targeted_elements_df):
    targeted_elements_df["grna_id_list"].str.split(",").explode("grna_id_list")
    return


@app.cell
def _():
    #exploded_grna_to_element = grna_to_element.assign(
    #    name=grna_to_element['name'].str.split(',')
    #).explode('name')

    #final_mapping = exploded_grna_to_element[['name', 'grna_target']].drop_duplicates().rename(columns={'name': 'grna_id'})

    # Use accession_id instead of grna_id
    #final_mapping = final_mapping.merge(
    #    guide_accession_map,
    #    left_on='grna_id',
    #    right_on='guide_id',
    #    how='left'
    #).drop(columns=['grna_id', 'guide_id']).rename(columns={'accession_id': 'grna_id'})[['grna_id', 'grna_target']]
    # write to a TSV file
    #final_mapping.to_csv(
    #    "results/03-2025-10-24_create_sceptre_input/grna_to_element_mapping_hg38.tsv",
    #    sep='\t',
    #    index=False,
    #    header=True,
    #    quoting=csv.QUOTE_NONE
    #)
    #final_mapping
    return


@app.cell
def _(mo):
    mo.md(
        r"""
    ## Create the discovery and positive pairs
    Now that we have the list of alignments we can see which ones fall in the regulatory region of a readout gene.
    """
    )
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
def _(genome_file, pybedtools, regulatory_region_pybed, targeted_elements_df):
    # Find all elements in the regulatory region for every gene (not just readout genes)
    # Require 50% of the element to overlap
    e_g_pairs = (
        regulatory_region_pybed
            .intersect(
                pybedtools.BedTool.from_dataframe(targeted_elements_df),
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
            .to_dataframe(header=None, usecols=[6,12])
            .rename(
                columns={6: 'response_id',
                         12: 'grna_target',
                         }
            )
    )

    discovery_pairs = (
        discovery_pairs
            .drop_duplicates()
            .reset_index(drop=True)
    )[["grna_target", "response_id"]]
    discovery_pairs.head()
    return (discovery_pairs,)


@app.cell
def _(
    e_g_pairs,
    genome_file,
    promoter_region_pybed,
    pybedtools,
    targeted_elements_df,
):
    # Find all elements in the regulatory region for every gene (not just readout genes)
    # Require 50% of the element to overlap
    e_tss_pairs = (
        promoter_region_pybed
            .intersect(
                pybedtools.BedTool.from_dataframe(targeted_elements_df),
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
            .to_dataframe(header=None, usecols=[6,12])
            .rename(
                columns={
                    12: 'grna_target',
                    6: 'response_id'
                         }
            )
    )
    # 
    positive_pairs = (
        positive_pairs
            .drop_duplicates()
            .reset_index(drop=True)
    )[["grna_target", "response_id"]]
    positive_pairs.head()
    return (positive_pairs,)


@app.cell
def _(mo):
    mo.md(r"""We exclude from the discovery pairs list the elements that overlap with TSS. They will be included as positive controls.""")
    return


@app.cell
def _(discovery_pairs, positive_pairs):
    # Remove all the pairs that overlaps with positive pairs
    keys = ['grna_target', 'response_id']

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
                readout_genes.rename(columns={"gene_id":"response_id"})[['response_id']],
                on='response_id'
            )[keys]
    )

    positive_pairs_final = (
        positive_pairs
            .merge(
                readout_genes.rename(columns={"gene_id":"response_id"})[['response_id']],
                on='response_id'
            )[keys]
    )

    discovery_pairs_final.to_csv(
        "results/03-2025-10-24_create_sceptre_input/discovery_pairs_no_pos_hg38_readout_genes.tsv",
        sep='\t',
        index=False,
        header=False,
        quoting=csv.QUOTE_NONE
    )

    positive_pairs_final.to_csv(
        "results/03-2025-10-24_create_sceptre_input/positive_pairs_hg38_readout_genes.tsv",
        sep='\t',
        index=False,
        header=False,
        quoting=csv.QUOTE_NONE
    )
    return discovery_pairs_final, positive_pairs_final


@app.cell
def _(discovery_pairs_final):

    discovery_pairs_final[["grna_target"]].value_counts()
    return


@app.cell
def _(positive_pairs_final):
    positive_pairs_final[["grna_target"]].value_counts()
    return


@app.cell
def _(discovery_pairs_final, positive_pairs_final):
    used_grna_targets = (
            set(
                discovery_pairs_final["grna_target"]
            ).union(
                positive_pairs_final["grna_target"]
            )
        )
    len(used_grna_targets)
    return (used_grna_targets,)


@app.cell
def _(mapped_guides_df):
    mapped_guides_df
    return


@app.cell
def _(csv, targeted_elements_df, used_grna_targets):
    targeting_guides = targeted_elements_df[
        targeted_elements_df['grna_target'].isin(used_grna_targets)
    ]["grna_id_list"].str.split(',').explode().unique()

    targeting_df = targeted_elements_df[
        targeted_elements_df['grna_target'].isin(used_grna_targets)
    ]

    exploded_grna_to_element = targeting_df.assign(
        grna_id=targeting_df['grna_id_list'].str.split(',')
    ).explode('grna_id')

    exploded_grna_to_element[["grna_id", "grna_target"]].drop_duplicates().to_csv("results/03-2025-10-24_create_sceptre_input/grna_to_element_df.tsv",
        sep='\t',
        index=False,
        header=True,
        quoting=csv.QUOTE_NONE)
    return


@app.cell
def _(targeted_elements_df, used_grna_targets):
    not_used_grna_targets = set(targeted_elements_df["grna_target"]) - set(used_grna_targets)
    # Number of elements not contributing to a discovery of positive pair.
    len(not_used_grna_targets)
    return (not_used_grna_targets,)


@app.cell
def _(all_guides, csv, not_used_grna_targets, pd, targeted_elements_df):
    safe_targeting_guides = targeted_elements_df[
        targeted_elements_df['grna_target'].isin(not_used_grna_targets)
    ]["grna_id_list"].str.split(',').explode().unique()

    non_targeting_guide = all_guides[all_guides["reason"].isin(["unmapped","discarded"])]

    non_targeting_sceptre_df = pd.DataFrame(
        {'grna_id': sorted(set(safe_targeting_guides).union(set(non_targeting_guide["guide_id"])))}
    )
    non_targeting_sceptre_df['grna_target'] = 'non-targeting'



    non_targeting_sceptre_df.to_csv("results/03-2025-10-24_create_sceptre_input/non_targeting_grna.tsv",
        sep='\t',
        index=False,
        header=False,
        quoting=csv.QUOTE_NONE)
    return


@app.cell(disabled=True, hide_code=True)
def _(genome_file, targeted_elements, tss_bed):
    closest_genes = targeted_elements.closest(tss_bed, d=True, g=genome_file)
    closest_genes.head()
    return (closest_genes,)


@app.cell(disabled=True, hide_code=True)
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
