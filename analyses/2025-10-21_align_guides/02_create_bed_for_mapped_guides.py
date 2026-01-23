import marimo

__generated_with = "0.18.0"
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
    # Guides not aligning to hg38
    unmapped_guide_fnp = "../../results/02-2025-10-21_align_guides/unmapped.tsv"
    # Discarded alignments
    discarded_alignments_fnp = "../../results/02-2025-10-21_align_guides/discarded.tsv"
    # Not valid alignments
    safe_targeting_fnp = "../../results/02-2025-10-21_align_guides/putative_safe_targeting.bed"
    # Mapped guides
    mapped_unique_guide_fnp = "../../results/02-2025-10-21_align_guides/guides_uniquely_mapped.bed"
    mapped_multi_guide_fnp = "../../results/02-2025-10-21_align_guides/guides_multimapping.bed"
    # Guides that are mapped after aligning to hg38 and intersecting with lifted over regions
    mapped_guide_overlapping_fnp = "../../results/02-2025-10-21_align_guides/guides_uniquely_mapped_overlapping_lifted_elements.bed"
    # Guides that are mapped after aligning to hg38 and intersecting with promoters and positive controls 
    mapped_overlapping_extra_fnp = "../../results/02-2025-10-21_align_guides/guides_uniquely_mapped_overlapping_list_of_promoters_and_positive_controls.bed"
    # Remaining mapped guides aligning to loci
    mapped_guides_custom_elements_alignments_fnp = "../../results/02-2025-10-21_align_guides/remaining_guides_in_loci.bed"
    # Remaining mapped guides custom elements
    mapped_guides_custom_elements_fnp = "../../results/02-2025-10-21_align_guides/remaining_guides_in_loci_merged_to_form_new_elements.bed"

    # Chromosome size file
    genome_file = "../../../../Annotations/HG38/GRCh38_EBV.chrom.sizes.no.alt.tsv"
    # TSS annotation file from GENCODE v44
    tss_file = "../../annotations/gencode_tss.bed"
    # Guide metadata file
    cellranger_feature_reference_file = "../../metadata/hiPSC-EC_feature-reference.csv"
    # DEPRECATED: Raw metadata guides. File unusable.
    #legacy_metadata_file = "metadata/240918_EC_TAPseq_CRISPR.full.tsv"
    # Gene metadata file
    gene_metadata_file = "../../metadata/metadata_genes_final.csv"
    return (
        cellranger_feature_reference_file,
        discarded_alignments_fnp,
        genome_file,
        mapped_guide_overlapping_fnp,
        mapped_guides_custom_elements_alignments_fnp,
        mapped_guides_custom_elements_fnp,
        mapped_overlapping_extra_fnp,
        mapped_unique_guide_fnp,
        safe_targeting_fnp,
        unmapped_guide_fnp,
    )


@app.cell
def _(mo):
    mo.md(r"""
    ## Load guide sequences metadata
    I am going to use the feature reference file passed as input to cellranger to map the old guide IDs to sequences. The legacy metadata has inconsistent accessions and contains redundancy and multiple accessions for the same sequence.
    I am assigning new accessions that are spacer based. Same spacer will give the same accession.
    """)
    return


@app.cell
def _(cellranger_feature_reference_file, csv, mo, pd):
    cellranger_feature_reference = pd.read_csv(cellranger_feature_reference_file, sep=',')
    cellranger_feature_reference['sequence'] = cellranger_feature_reference['sequence'].str.upper().str[1:]

    message = (
        f"Loaded {len(cellranger_feature_reference)} guide sequences "
        "from CellRanger feature reference."
    )
    mo.md(message)
    guide_metadata_df = cellranger_feature_reference[['sequence', 'id']].rename(columns={'sequence':'guide_id', 'id': 'aliases' })

    # Save the final guide metadata with consistent accessions
    guide_metadata_df.to_csv(
        "../../results/03-2025-10-24_create_sceptre_input/guide_metadata_with_new_accessions.tsv",
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
    mo.md(r"""
    ## Load the alignment results
    After re-aligning and annotating manually the guides to hg38 I am loading all the different categories of alignments to create the final list of discovery pairs.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ### Unmapped guides
    These are the guides that did not align to hg38 at all.
    """)
    return


@app.cell
def _(alias_to_guide_id, pd, unmapped_guide_fnp):
    # Load the unmapped guides
    unmapped_guides = pd.read_csv(
        unmapped_guide_fnp,
        sep='\t',
        header=0
    )
    unmapped_guides['alias'] = unmapped_guides['read_name']
    unmapped_guides['guide_id'] = unmapped_guides['read_name'].map(alias_to_guide_id)
    unmapped_guides.head()
    return (unmapped_guides,)


@app.cell
def _(mo):
    mo.md(r"""
    ### Putative safe targeting guides
    These are the guides that aligned to hg38 but do not overlap with any of the 16 designed loci.
    We need to pay attention to make sure that thtese guides are not overlapping promoters or positive controls.
    """)
    return


@app.cell
def _(alias_to_guide_id, pd, safe_targeting_fnp):
    # Loading the putative safe targeting guides
    safe_targeting_bed = pd.read_csv(safe_targeting_fnp, sep='\t', header=None, names=['chrom', 'start', 'end', 'guide_id', 'score', 'strand', 'NM', 'AS'])
    safe_targeting_bed['alias'] = safe_targeting_bed['guide_id']
    safe_targeting_bed['guide_id'] = safe_targeting_bed['guide_id'].map(alias_to_guide_id)
    safe_targeting_bed.head()
    return


@app.cell
def _(mo):
    mo.md(r"""
    ### Discarded alignments
    These are the guides that aligned to hg38 but were discarded because the PAM wasn't aligned properly or there were too many mismatches.
    Some of the guides in this category might have valid alignments in other regions. So, to find how many guides were considered unmapped we need to remove the guides in the mapped category that are also in the discarded category.
    """)
    return


@app.cell
def _(alias_to_guide_id, discarded_alignments_fnp, pd):
    # Loading the discarded alignments
    discarded_alignments = pd.read_csv(
        discarded_alignments_fnp,
        sep='\t',
        header=0
    )
    discarded_alignments['alias'] = discarded_alignments['read_name']
    discarded_alignments['guide_id'] = discarded_alignments['read_name'].map(alias_to_guide_id)
    discarded_alignments.head()
    return (discarded_alignments,)


@app.cell
def _(mo):
    mo.md(r"""
    ### Mapped guides
    The guides that aligned to hg38 can be split into multiple subcategories:

    `mapped_unique_guide_fnp`: This is the complete set of uniquely mapped guides. 'Uniquely mapped' indicates that either the guides align to a single location in the genome or that other alignments were discarded due to poor PAM or too many mismatches.

    `mapped_guide_overlapping_fnp`: This is a subset of the uniquely mapped guides that overlap the designed elements after lifting over.

    `mapped_overlapping_extra_fnp`: This is the set of guide that do not overlap the designed elements but overlaps promoters and positive controls from previous experiments.

    `mapped_guides_custom_elements_alignments_fnp`: These are the guides that aligned within the boundaries of the 16 designed loci but do not overlap any lifted elements, promoters or positive controls. This is most likely due to small differences in the coordinates after liftover. These guides will be further processed to merge their alignments and form new custom elements within the loci.

    `mapped_guides_custom_elements_fnp`. These are the new custom elements formed by merging the alignments of the remaining guides in loci.

    `mapped_multi_guide_fnp`. Guides with more than one valid alignment. These guides will be discarded from the final list.
    """)
    return


@app.cell
def _(alias_to_guide_id, pd):
    # Load all the mapped guides
    mapped_bed = pd.read_csv(
        "../../results/02-2025-10-21_align_guides/mapped.bed",
        sep='\t',
        header=None,
        names=['chrom', 'start', 'end', 'guide_id', 'score', 'strand', 'NM', 'AS']
    )
    mapped_bed['alias'] = mapped_bed['guide_id']
    mapped_bed['guide_id'] = mapped_bed['guide_id'].map(alias_to_guide_id)
    mapped_bed.head()
    return (mapped_bed,)


@app.cell
def _(
    alias_to_guide_id,
    discarded_alignments,
    mapped_bed,
    pd,
    unmapped_guides,
):
    # Creating the final list of unmapped guides
    # Guides without a single valid alignment
    no_valid_alignments = set(discarded_alignments.alias) - set(mapped_bed.alias)
    final_unmapped_guides = pd.DataFrame(no_valid_alignments, columns=['alias'])
    final_unmapped_guides['guide_id'] = final_unmapped_guides['alias'].map(alias_to_guide_id)
    final_unmapped_guides = final_unmapped_guides[['guide_id', 'alias']]
    final_unmapped_guides['reason'] = 'no_valid_alignment'
    # Add the unmapped guides
    final_unmapped_guides = pd.concat([
        final_unmapped_guides,
        unmapped_guides[['guide_id', 'alias']].assign(reason='no_alignment')
    ], ignore_index=True)
    # Write to file
    final_unmapped_guides.to_csv(
        "../../results/02-2025-10-21_align_guides/non_targeting_guides_final.tsv",
        sep='\t',   
        index=False
    )
    final_unmapped_guides.head()
    return


@app.cell
def _(alias_to_guide_id, mapped_unique_guide_fnp, pd):
    # Loading the uniquely mapped guides
    mapped_unique_bed = pd.read_csv(
        mapped_unique_guide_fnp,
        sep='\t',
        header=None,
        names=['chrom', 'start', 'end', 'guide_id', 'score', 'strand', 'NM', 'AS']
    )
    mapped_unique_bed['alias'] = mapped_unique_bed['guide_id']
    mapped_unique_bed['guide_id'] = mapped_unique_bed['guide_id'].map(alias_to_guide_id)

    # Write to file
    mapped_unique_bed.to_csv(
        "../../results/02-2025-10-21_align_guides/uniquely_mapped_guides_final.bed",
        sep='\t',   
        index=False
    )
    mapped_unique_bed.head()
    return


@app.cell
def _(mo):
    mo.md(r"""
    ### OUTPUTS
    The relevant outputs from this notebook are:
    - `uniquely_mapped_guides_final.bed`: This is the final list of guides that are uniquely mapped to hg38 (n=43235).
    - `non_targeting_guides_final.tsv`: This is the list of guides that either don't map to hg38 or don't have a valid mapping (n= 501).

    The remaining guide to reach the 43778 total are the ones that map to multiple locations and will be discarded (n=42).
    """)
    return


@app.cell
def _(alias_to_guide_id, mapped_guide_overlapping_fnp, pd):
    # Loading the guides that mapped to lifted elements
    mapped_guide_overlapping_bed = pd.read_csv(mapped_guide_overlapping_fnp, sep='\t', header=None,
    names=['chrom_element', 'start_element', 'end_element', 'chrom', 'start', 'end', 'guide_id', 'score', 'strand', 'NM', 'AS'])
    mapped_guide_overlapping_bed['alias'] = mapped_guide_overlapping_bed['guide_id']
    mapped_guide_overlapping_bed['guide_id'] = mapped_guide_overlapping_bed['guide_id'].map(alias_to_guide_id)
    mapped_guide_overlapping_bed.head()
    return


@app.cell
def _(alias_to_guide_id, mapped_overlapping_extra_fnp, pd):
    # Loading the guides mapping to TSS and positive controls
    mapped_guide_overlapping_extra_bed = pd.read_csv(mapped_overlapping_extra_fnp, sep='\t', header=None,
    names=['chrom', 'start', 'end', 'guide_id', 'score', 'strand', 'NM', 'AS', 'chrom_element', 'start_element', 'end_element', 'name_element', 'alias_element', 'pool_element'])
    mapped_guide_overlapping_extra_bed['alias'] = mapped_guide_overlapping_extra_bed['guide_id']  
    mapped_guide_overlapping_extra_bed['guide_id'] = mapped_guide_overlapping_extra_bed['guide_id'].map(alias_to_guide_id)
    mapped_guide_overlapping_extra_bed.head()
    return


@app.cell
def _(
    alias_to_guide_id,
    mapped_guides_custom_elements_alignments_fnp,
    mapped_guides_custom_elements_fnp,
    pd,
):
    # Guides aligned to remaining custom elements
    mapped_guides_custom_elements_alignments_bed = pd.read_csv(mapped_guides_custom_elements_alignments_fnp, sep='\t', header=None,
    names=['chrom', 'start', 'end', 'guide_id', 'score', 'strand', 'NM', 'AS', 'chrom_locus', 'start_locus', 'end_locus', 'central_gene']
    )
    mapped_guides_custom_elements_alignments_bed['alias'] = mapped_guides_custom_elements_alignments_bed['guide_id']
    mapped_guides_custom_elements_alignments_bed['guide_id'] = mapped_guides_custom_elements_alignments_bed['guide_id'].map(alias_to_guide_id)
    mapped_guides_custom_elements_alignments_bed.head()

    # Load the new elements with the guides mapping to them
    new_elements_bed = pd.read_csv(mapped_guides_custom_elements_fnp, sep='\t', header=None,
    names=['chrom', 'start', 'end', 'guide_mapping_list']
    )

    new_elements_bed.head()
    return


@app.cell
def _():
    return


@app.cell
def _(discarded_guides, pd, unmapped_guides):
    all_guides = pd.concat([
        pd.DataFrame(list(discarded_guides), columns=['guide_id']).assign(reason='discarded'),
        unmapped_guides
    ])
    all_guides.groupby('reason')['guide_id'].count()
    return (all_guides,)


@app.cell
def _(all_guides, csv):
    all_guides.to_csv("results/03-2025-10-24_create_sceptre_input/non_targeting_grna.tsv",
        sep='\t',
        index=False,
        header=True,
        quoting=csv.QUOTE_NONE)
    return


@app.cell
def _(mo):
    mo.md(r"""
    Elements from lifted over design file
    """)
    return


@app.cell
def _(csv, valid_overlapping):
    designed_elements = valid_overlapping[['element_chrom', 'element_start', 'element_end', 'element_name', 'alias']].drop_duplicates().reset_index(drop=True)
    designed_elements_df = designed_elements.rename(columns={
        'element_chrom': 'chrom',
        'element_start': 'start',
        'element_end': 'end',
        'element_name': 'designed_element_name'
    })
    designed_elements_df['designed_element_size'] = designed_elements_df['end'] - designed_elements_df['start']
    designed_elements_df.to_csv("results/03-2025-10-24_create_sceptre_input/designed_elements.tsv",
            sep="\t",
            index=False,
            header=True,
            quoting=csv.QUOTE_NONE,
                               )
    designed_elements_df
    return (designed_elements_df,)


@app.cell
def _(designed_elements_df):
    # Number of elements
    len(designed_elements_df)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## Creating the discovery and positive pairs
     - Intersect the valid mapping with the regulatory regions of the readout genes.
     - Intersect the valid mapping guides with the readout genes TSSs.
     - Intersect the safe_targeting with the regulatory regions of the readout genes and exclude the ones that overlap.
    """)
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
def _(designed_elements_df, genome_file, pybedtools, regulatory_region_pybed):
    # Find all elements in the regulatory region for every gene (not just readout genes)
    # Require 50% of the element to overlap
    e_g_pairs = (
        regulatory_region_pybed
            .intersect(
                pybedtools.BedTool.from_dataframe(designed_elements_df),
                g=genome_file,
                wa=True,
                wb=True,
                F=0.5,
            )
    )
    e_g_pairs.head()
    return (e_g_pairs,)


@app.cell
def _(readout_genes):
    # Filter discovery_pairs and keep only readout genes
    readout_gene_ids = set(readout_genes['gene_id'])
    readout_gene_ids
    return (readout_gene_ids,)


@app.cell
def _(e_g_pairs, readout_gene_ids):
    # Format the results in a dataframe
    discovery_pairs = (
        e_g_pairs
            .to_dataframe(header=None, usecols=[6, 7, 8, 12, 13])
            .rename(
                columns={6: 'response_id',
                         7: 'gene_type',
                         8: 'gene_symbol',
                         12: 'grna_target',
                         13: 'grna_target_alias'
                         }
            )
    )

    discovery_pairs = (
        discovery_pairs
            .drop_duplicates()
            .reset_index(drop=True)
    )[["grna_target", "response_id", "gene_symbol", "gene_type", "grna_target_alias"]] # Reordering columns.
    discovery_pairs = discovery_pairs[
        discovery_pairs['response_id'].isin(readout_gene_ids)
    ].reset_index(drop=True)
    discovery_pairs.head(), len(discovery_pairs)
    return (discovery_pairs,)


@app.cell
def _(
    designed_elements_df,
    e_g_pairs,
    genome_file,
    promoter_region_pybed,
    pybedtools,
):
    # Find all elements in the regulatory region for every gene (not just readout genes)
    # Require 50% of the element to overlap
    e_tss_pairs = (
        promoter_region_pybed
            .intersect(
                pybedtools.BedTool.from_dataframe(designed_elements_df),
                g=genome_file,
                wa=True,
                wb=True,
                F=0.5,
            )
    )
    e_g_pairs.head()
    return (e_tss_pairs,)


@app.cell
def _(e_tss_pairs, readout_gene_ids):
    # Format the results in a dataframe
    positive_pairs = (
        e_tss_pairs
            .to_dataframe(header=None, usecols=[6, 7, 8, 12, 13])
            .rename(
                columns={6: 'response_id',
                         7: 'gene_type',
                         8: 'gene_symbol',
                         12: 'grna_target',
                         13: 'grna_target_alias'
                         }
            )
    )
    # 
    positive_pairs = (
        positive_pairs
            .drop_duplicates()
            .reset_index(drop=True)
    )[["grna_target", "response_id", "gene_symbol", "gene_type", "grna_target_alias"]]
    positive_pairs = positive_pairs[
        positive_pairs['response_id'].isin(readout_gene_ids)
    ].reset_index(drop=True)
    positive_pairs, len(positive_pairs)
    return (positive_pairs,)


@app.cell
def _(
    csv,
    discovery_pairs_final,
    gene_id_name_map,
    pd,
    positive_pairs_final,
    readout_gene_ids,
):
    missing_genes = set(readout_gene_ids) - set(discovery_pairs_final["response_id"]).union(set(positive_pairs_final["response_id"]))
    # Write the list to file adding gene symbols
    missing_genes_df = pd.DataFrame({'gene_id': sorted(missing_genes)})
    missing_genes_df = (
        missing_genes_df
            .merge(
                gene_id_name_map,
                on='gene_id',
                how='left'
            )
    )
    missing_genes_df.to_csv("results/03-2025-10-24_create_sceptre_input/missing_readout_genes.tsv",
        sep='\t',
        index=False,
        header=True,
        quoting=csv.QUOTE_NONE)
    missing_genes_df.head()
    return


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
def _(mo):
    mo.md(r"""
    We exclude from the discovery pairs list the elements that overlap with TSS. They will be included as positive controls.
    """)
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
    return


@app.cell
def _(csv, discovery_pairs, positive_pairs):
    discovery_pairs.to_csv(
        "results/03-2025-10-24_create_sceptre_input/discovery_pairs_no_pos_hg38_readout_genes.tsv",
        sep='\t',
        index=False,
        header=False,
        quoting=csv.QUOTE_NONE
    )

    positive_pairs.to_csv(
        "results/03-2025-10-24_create_sceptre_input/positive_pairs_hg38_readout_genes.tsv",
        sep='\t',
        index=False,
        header=False,
        quoting=csv.QUOTE_NONE
    )
    return


@app.cell
def _(discovery_pairs):

    discovery_pairs[["grna_target"]].value_counts()
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


@app.cell(disabled=True)
def _(mo):
    mo.md(r"""
    Load the BED file with the results of the alignment. Convert the names to accessions and create the pybedtools object.
    """)
    return


@app.cell(disabled=True)
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


@app.cell(disabled=True, hide_code=True)
def _():
    #def _():
    #    all_guides = pd.concat([
    #        pd.DataFrame(list(kept_guides), columns=#['guide_id']).assign(reason='kept'),
    #        pd.DataFrame(list(discarded_guides), columns=['guide_id']).assign(reason='discarded'),
    #        unmapped_guides
    #    ])
    #    return all_guides.groupby('reason')['guide_id'].count()


    #_()
    return


@app.cell(disabled=True)
def _(mo):
    mo.md(r"""
    ## Create list of assayed elements
    Now that we have the coordinates of the guides, we can merge the ones close to each other to create the list of elements being targeted by the same guides.
    Parameters used in the analysis:

     - guide_merge_threshold = 50 # Merge together guides closer than dist_threshold
     - promoter_distance = 500 # Define the promoter as +/- 500bp from TSS
     - regulatory_region_distance = 1_000_000 # Extend +/- 1Mb around the TSS
    """)
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


if __name__ == "__main__":
    app.run()
