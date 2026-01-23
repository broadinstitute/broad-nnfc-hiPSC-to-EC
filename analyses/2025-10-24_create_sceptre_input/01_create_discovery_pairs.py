import marimo

__generated_with = "0.18.1"
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
    return csv, pd, pybedtools


@app.cell
def _():
    promoter_distance = 500 # Define the promoter as +/- 500bp from TSS
    regulatory_region_distance = 1_000_000 # Extend +/- 1Mb around the TSS
    return promoter_distance, regulatory_region_distance


@app.cell
def _():
    # Loading the element coordinates I want to annotate
    elements_fnp = "../../results/03-2025-10-24_create_sceptre_input/all_element_annotations_final.bed"

    # Chromosome size file
    genome_file = "../../../../Annotations/HG38/GRCh38_EBV.chrom.sizes.no.alt.tsv"
    # TSS annotation file from GENCODE v44
    tss_file = "../../annotations/gencode_tss.bed"
    # Guide metadata file
    cellranger_feature_reference_file = "../../metadata/hiPSC-EC_feature-reference.csv"
    # Genes metadata file
    metadata_genes_fnp = "../../metadata/metadata_genes_final.csv"
    return elements_fnp, genome_file, metadata_genes_fnp, tss_file


@app.cell
def _(mo):
    mo.md(r"""
    ## Load elements
    Loading the designed elements file and converting it to a pybedtools object.
    """)
    return


@app.cell
def _(elements_fnp, pd, pybedtools):
    elements_df = pd.read_csv(
        elements_fnp,
        sep="\t",
        header=None,
        names=[
            "chrom",
            "start",
            "end",
            "element_id",
            "source",
        ],
    )
    elements_pybed = pybedtools.BedTool.from_dataframe(elements_df)
    elements_df.head(), elements_pybed.to_dataframe().head()
    return (elements_pybed,)


@app.cell
def _(mo):
    mo.md(r"""
    ## Load genes metadata
    Loading the list of readout genes and the TSS for all the genes in GENCODE v44.
    """)
    return


@app.cell
def _(metadata_genes_fnp, pd):
    # List of promoters with validated guides
    promoters_with_validated_guides = ['CDH5', 'CTNND1', 'DLL4', 'FAM136A', 'FJX1', 'HAND1', 'HEY2', 'KDR', 'MOGS', 'PRKD1', 'PYCARD', 'RBMXL1', 'SAT2', 'SOX2', 'SRPK1', 'TMEM107', 'TMEM256', 'VKORC1']

    # Get the list of readout genes
    readout_genes = pd.read_csv(
        metadata_genes_fnp
    )
    readout_genes =(
        readout_genes["intended_target_gene_symbol"]
        .drop_duplicates()
        .reset_index(drop=True)
    )
    readout_genes = readout_genes.to_frame()
    readout_genes.head(), readout_genes.shape
    return promoters_with_validated_guides, readout_genes


@app.cell
def _(genome_file, pybedtools, tss_file):
    # Loading the gene TSS annotation
    tss_bed = pybedtools.BedTool(tss_file).sort(g=genome_file)
    tss_bed.head()
    return (tss_bed,)


@app.cell
def _(promoters_with_validated_guides, readout_genes):
    # Check if all the genes with validated promoters are in the readout genes
    set(readout_genes["intended_target_gene_symbol"]).intersection(promoters_with_validated_guides)
    missing_positive_control_genes = set(promoters_with_validated_guides) - set(readout_genes["intended_target_gene_symbol"])
    if len(missing_positive_control_genes) > 0:
        print(f"Warning: The following genes with validated promoters are missing from the readout genes: {missing_positive_control_genes}")
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## Creating the discovery and positive pairs
     - For each readout gene, find all the elements within +/- 1Mb of its TSS.
     - For each readout gene, find all the elements overlapping its TSS +/- 500bp (positive control pairs).
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
def _(mo):
    mo.md(r"""
    ### Discovery pairs
    """)
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
def _(elements_pybed, genome_file, regulatory_region_pybed):
    # Find all elements in the regulatory region for every gene (not just readout genes)
    # Require 50% of the element to overlap
    discovery_pairs_all = (
        regulatory_region_pybed
            .intersect(
                elements_pybed,
                g=genome_file,
                wa=True,
                wb=True,
                F=0.5,
            )
    )
    discovery_pairs_all.head()
    # Write out the unfilterd discovery pairs
    discovery_pairs_all_df = discovery_pairs_all.to_dataframe(names=[
        "regulatory_region_chrom",
        "regulatory_region_start",
        "regulatory_region_end",
        "transcript_id",
        "score",
        "strand",
        "gene_id",
        "gene_type",
        "gene_symbol",
        "element_chrom",
        "element_start",
        "element_end",
        "element_id",
        "element_source",
    ])
    discovery_pairs_all_df.to_csv(
        "../../results/03-2025-10-24_create_sceptre_input/discovery_pairs_all.bed",
        sep="\t",
        header=True,
        index=False,
    )
    return (discovery_pairs_all_df,)


@app.cell
def _(discovery_pairs_all_df, readout_genes):
    # Keep only discovery pairs for readout genes
    discovery_pairs_readout_genes = discovery_pairs_all_df[
        discovery_pairs_all_df["gene_id"].isin(readout_genes["gene_id"])
    ]
    # Write out the filtered discovery pairs
    discovery_pairs_readout_genes.to_csv(
        "../../results/03-2025-10-24_create_sceptre_input/discovery_pairs_readout_genes.bed",
        sep="\t",
        header=True,
        index=False,
    )
    discovery_pairs_readout_genes.head()
    return (discovery_pairs_readout_genes,)


@app.cell
def _(discovery_pairs_readout_genes, pd, readout_genes):
    # Find which readout genes are missing from discovery pairs
    missing_pairs_for_readout_genes = (
        set(readout_genes.intended_target_gene_symbol)
        - set(discovery_pairs_readout_genes.gene_symbol)
    )

    # Convert set → Series
    missing_pairs_for_readout_genes = pd.Series(sorted(missing_pairs_for_readout_genes))
    missing_pairs_for_readout_genes.name = "readout_genes_without_discovery_pairs"

    # Save to file
    missing_pairs_for_readout_genes.to_csv(
        "../../results/03-2025-10-24_create_sceptre_input/missing_discovery_pairs_for_readout_genes.txt",
        sep="\t",
        header=True,
        index=False,
    )

    missing_pairs_for_readout_genes

    return


@app.cell
def _(mo):
    mo.md(r"""
    ### Promoter pairs (positive controls)
    """)
    return


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
def _(elements_pybed, genome_file, promoter_region_pybed):
    # Find all elements in the promoter region for every gene (not just readout genes)
    # Require 50% of the element to overlap
    promoter_pairs_readout_genes_bed = (
        promoter_region_pybed
            .intersect(
                elements_pybed,
                g=genome_file,
                wa=True,
                wb=True,
                F=0.5,
            )
    )
    promoter_pairs_readout_genes = promoter_pairs_readout_genes_bed.to_dataframe(names=[
        "promoter_region_chrom",
        "promoter_region_start",
        "promoter_region_end",
        "transcript_id",
        "score",
        "strand",
        "gene_id",
        "gene_type",
        "gene_symbol",
        "element_chrom",
        "element_start",
        "element_end",
        "element_id",
        "element_source",
    ])
    # Save promoter-element pairs to file
    promoter_pairs_readout_genes.to_csv(
        "../../results/03-2025-10-24_create_sceptre_input/promoter_element_pairs_all.bed",
        sep="\t",
        header=True,
        index=False,
    )

    promoter_pairs_readout_genes.head()
    return (promoter_pairs_readout_genes,)


@app.cell
def _(pd, promoter_pairs_readout_genes, readout_genes):
    # Find which readout genes are missing from promoter pairs
    missing_tss_pairs_for_readout_genes = (
        set(readout_genes.intended_target_gene_symbol)
        - set(promoter_pairs_readout_genes.gene_symbol)
    )

    # Convert set → Series
    missing_tss_pairs_for_readout_genes = pd.Series(sorted(missing_tss_pairs_for_readout_genes))
    missing_tss_pairs_for_readout_genes.name = "readout_genes_without_promoter_pairs"

    # Save to file
    missing_tss_pairs_for_readout_genes.to_csv(
        "../../results/03-2025-10-24_create_sceptre_input/missing_promoter_pairs_for_readout_genes.txt",
        sep="\t",
        header=True,
        index=False,
    )

    missing_tss_pairs_for_readout_genes


    return


@app.cell
def _(
    csv,
    discovery_pairs_readout_genes,
    gene_id_name_map,
    pd,
    promoter_pairs_readout_genes,
    readout_genes,
):
    missing_genes = set(readout_genes["gene_id"]) - set(discovery_pairs_readout_genes["gene_id"]).union(set(promoter_pairs_readout_genes["gene_id"]))
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
    missing_genes_df.to_csv("../../results/03-2025-10-24_create_sceptre_input/missing_readout_genes.tsv",
        sep='\t',
        index=False,
        header=True,
        quoting=csv.QUOTE_NONE)
    missing_genes_df.head(), missing_genes_df.shape
    return


@app.cell
def _(mo):
    mo.md(r"""
    We exclude from the discovery pairs list the elements that overlap with TSS. They will be included as positive controls.
    """)
    return


@app.cell
def _(discovery_pairs_readout_genes, promoter_pairs_readout_genes):
    # Remove all the pairs that overlaps with positive pairs

    keys = ['element_id', 'gene_id', 'gene_symbol']
    dedup_promoter_pairs_readout_genes = promoter_pairs_readout_genes[keys].drop_duplicates(keep='first')

    promoter_pairs_readout_genes[[*keys,"element_source"]].drop_duplicates(keep='first').to_csv(
        "../../results/03-2025-10-24_create_sceptre_input/positive_control_pairs_dedup_final.tsv",
        sep="\t",
        header=["grna_target", "response_id","gene_symbol", "element_source"],
        index=False,
    )

    dedup_discovery_pairs_readout_genes = discovery_pairs_readout_genes[keys].drop_duplicates(keep='first')
    dedup_discovery_pairs_readout_genes_no_pos = (
        dedup_discovery_pairs_readout_genes
            .merge(
                dedup_promoter_pairs_readout_genes,
                on=keys,
                how='left',
                indicator=True
            )
            .query('_merge == "left_only"')
            .drop(columns=['_merge'])
    )
    dedup_discovery_pairs_readout_genes_no_pos.to_csv(
        "../../results/03-2025-10-24_create_sceptre_input/discovery_pairs_no_positive_controls_final.tsv",
        sep="\t",
        header=["grna_target", "response_id","gene_symbol"],
        index=False,
    )

    dedup_discovery_pairs_readout_genes_no_pos.head(), dedup_discovery_pairs_readout_genes_no_pos.shape
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## Final files
    `all_guides_annotations_final.tsv`: File containing the mapping of all guides to elements.
    `discovery_pairs_no_positive_controls_final.tsv`: File containing all discovery pairs (element, readout gene) excluding positive controls.
    `positive_control_pairs_dedup_final.tsv`: File containing all positive control pairs (element, readout gene).
    """)
    return


if __name__ == "__main__":
    app.run()
