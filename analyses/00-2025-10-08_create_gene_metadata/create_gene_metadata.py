import marimo

__generated_with = "0.15.2"
app = marimo.App()


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell
def _(mo):
    mo.md(r"This notebook reads the results of the primer alignment to transcript FASTA file and check that the gene we associated with the primer is indeed present in the transcript FASTA file.")
    return


@app.cell
def _():
    import pandas as pd
    return (pd,)


@app.cell
def _(mo):
    mo.md(r"Load the raw gene metadata file.")
    return


@app.cell
def _(pd):
    metadata_raw = pd.read_csv("metadata/metadata_genes_raw.csv")
    metadata_raw
    return (metadata_raw,)


@app.cell
def _(mo):
    mo.md(r"Load the primer mapping results.")
    return


@app.cell
def _(pd):
    map_results = pd.read_csv("results/2025-10-08_create_gene_metadata/intermediates/primer_matches_results.csv")
    map_results

    return (map_results,)


@app.cell
def _(mo):
    mo.md(r"Create a gene symbol to gene ID mapping using the results of the primer mapping.")
    return


@app.cell
def _(map_results):
    gene_symbol_to_id = {row["gene_symbol"]: row["gene_id"] for _, row in map_results.iterrows()}
    gene_symbol_to_id
    return (gene_symbol_to_id,)


@app.cell
def _(mo):
    mo.md(r"List the intended target genes")
    return


@app.cell
def _(map_results):
    intended_targets = map_results["intended_target_gene_symbol"].unique().tolist()
    intended_targets
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        Find all the primers for which the intended target gene is not present in the transcript FASTA file.
        Create a dictionary mapping the intended target gene to the actual gene symbol found in the transcript FASTA file.
        """
    )
    return


@app.cell
def _(map_results):
    # Display the rows with only non-matching intended targets
    match_counts = map_results.groupby("inner_primer")["is_intended_target"].sum()

    # Display the non-matching rows
    mask_not_matching = map_results["inner_primer"].isin(match_counts[match_counts == 0].index)

    gene_name_update_dict = {row["intended_target_gene_symbol"]: row["gene_symbol"] for _, row in map_results[mask_not_matching].iterrows()}
    gene_name_update_dict
    return (gene_name_update_dict,)


@app.cell
def _(mo):
    mo.md(r"Load the central and locus genes annotations. I am also updating the gene names if needed and save the new files.")
    return


@app.cell
def _(gene_name_update_dict, pd):
    locus_annotations = pd.read_csv("metadata/ec_screen_16_loci_2mb_info_raw.tsv", sep="\t", header=0)
    locus_bed = pd.read_csv("metadata/ec_screen_16_loci_2mb_locus_raw.bed", sep="\t", header=None, names=["chrom", "start", "end", "gene_symbol"])

    # Update the columns `central_gene`, `genes_above_tpm_ids`, `ubiq_genes`, and `gene_names` using the gene_name_update_dict
    locus_annotations["central_gene"] = locus_annotations["central_gene"].replace(gene_name_update_dict)
    locus_annotations["gene_names"] = locus_annotations["gene_names"].apply(    lambda x: ",".join([gene_name_update_dict.get(gene, gene) for gene in x.split(",")]) if pd.notna(x) else x
    )
    locus_annotations["genes_above_tpm_ids"] = locus_annotations["genes_above_tpm_ids"].apply(    lambda x: ",".join([gene_name_update_dict.get(gene, gene) for gene in x.split(",")]) if pd.notna(x) else x
    )
    locus_annotations["ubiq_genes"] = locus_annotations["ubiq_genes"].apply(    lambda x: ",".join([gene_name_update_dict.get(gene, gene) for gene in x.split(",")]) if pd.notna(x) else x
    )
    locus_bed["gene_symbol"] = locus_bed["gene_symbol"].replace(gene_name_update_dict)
    # Update column names
    locus_annotations = locus_annotations.rename(columns={
        "locus": "central_gene_id", 
        "central_gene": "central_gene_symbol",
        "gene_names": "locus_gene_symbols",
        "genes_above_tpm_ids": "gene_above_tpm_symbols",
        "n_ubiq_genes": "n_ubiquitinous_genes",
        "ubiq_genes": "ubiquitinous_gene_symbols"
    })
    locus_annotations.to_csv("metadata/ec_screen_16_loci_2mb_info_final.tsv", sep="\t", index=False)
    locus_bed.to_csv("metadata/ec_screen_16_loci_2mb_locus_final.bed", sep="\t", index=False, header=False)
    locus_annotations
    return (locus_annotations,)


@app.cell
def _(mo):
    mo.md(r"For each gene in `central_gene` and `gene_names` columns create a mapping to a tuple (is_central_gene: bool, locus_central_gene: gene_symbol).")
    return


@app.cell
def _(locus_annotations):
    gene_central_locus_mapping = {}
    total_genes_added = 0
    locus_annotations["central_gene_symbol"]
    for _, row in locus_annotations.iterrows():
        central_gene = row["central_gene_symbol"]
        gene_central_locus_mapping[central_gene] = (True, central_gene)
        total_genes_added += 1
        for locus_gene in row["locus_gene_symbols"].split(";"):
            if gene_central_locus_mapping.get(locus_gene, [False])[0]:
                continue
            gene_central_locus_mapping[locus_gene] = (False, central_gene)
            total_genes_added += 1


    total_genes_added,gene_central_locus_mapping
    return (gene_central_locus_mapping,)


@app.cell
def _(gene_central_locus_mapping):
    gene_central_locus_mapping["BMPR1A"]
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        Generate the new metadata files with the corrected gene symbols.
    
        ### metadata_genes_final.csv
    
        Columns:
    
         - intended_target_gene_id: Ensembl gene ID
         - intended_target_gene_symbol: Gene symbol
         - inner_primer: Inner primer sequence
         - outer_primer: Outer primer sequence
         - inner_primer_adapter: Inner primer sequence with adapters
         - is_central_gene: Is central gene (True/False)
         - locus_central_gene: Locus central gene the gene belongs to (if applicable)
         - well: Well position they came from the commercial provider
         - batch: Batch information
        """
    )
    return


@app.cell
def _(
    gene_central_locus_mapping,
    gene_name_update_dict,
    gene_symbol_to_id,
    metadata_raw,
):
    import datetime

    metadata_final = metadata_raw
    # Update the gene symbols
    metadata_final["intended_target_gene_symbol"] = metadata_final["intended_target_gene_symbol"].replace(gene_name_update_dict)
    # Update the gene IDs
    metadata_final["intended_target_gene_id"] = metadata_final["intended_target_gene_symbol"].replace(gene_symbol_to_id)
    # Add boolean column if the gene is a central gene and the locus central gene
    metadata_final["is_central_gene"] = metadata_final["intended_target_gene_symbol"].apply(
        lambda x: gene_central_locus_mapping[x][0] if x in gene_central_locus_mapping else False
    )
    # Add the locus central gene symbol
    metadata_final["locus_central_gene"] = metadata_final["intended_target_gene_symbol"].apply(
        lambda x: gene_central_locus_mapping[x][1] if x in gene_central_locus_mapping else None
    )

    # Reorder columns
    metadata_final = metadata_final[["intended_target_gene_id","intended_target_gene_symbol","outer_primer","inner_primer","inner_primer_adapter","is_central_gene",
        "locus_central_gene",
        "well",
        "batch"
    ]]



    # Save the final metadata
    metadata_final.to_csv("metadata/metadata_genes_final.csv", index=False)

    metadata_final
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
