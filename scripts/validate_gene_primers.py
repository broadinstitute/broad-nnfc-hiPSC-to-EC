import marimo

__generated_with = "0.15.2"
app = marimo.App()


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell
def _(mo):
    mo.md(r"""This notebook reads the results of the primer alignment to transcript FASTA file and check that the gene we associated with the primer is indeed present in the transcript FASTA file.""")
    return


@app.cell
def _():
    import pandas as pd
    return (pd,)


@app.cell
def _(pd):
    map_results = pd.read_csv("metadata/hiPSC-EC/primer_matches.csv")
    map_results
    return (map_results,)


@app.cell
def _(map_results, mo):
    # Group by inner primer and make  sure there is at least one 'is_intended_target' to true
    grouped = map_results.groupby("inner_primer")
    intended_targets = grouped["is_intended_target"].sum()
    intended_targets[intended_targets == 0]
    # If there are any inner primers with no intended targets, raise an error

    mo.md(r"""All inner primers have at least one intended target. The primer design is consistent with the transcript FASTA file.""")
    intended_targets[intended_targets == 0]
    return (intended_targets,)


@app.cell
def _(intended_targets, map_results):
    # Display the rows with no intended targets
    mask_missing_genes = map_results["inner_primer"].isin(intended_targets[intended_targets == 0].index)
    map_results[["inner_primer", "intended_target_gene_symbol", "gene_symbol"]][mask_missing_genes]
    return


if __name__ == "__main__":
    app.run()
