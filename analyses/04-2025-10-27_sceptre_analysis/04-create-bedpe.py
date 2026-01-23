import marimo

__generated_with = "0.18.1"
app = marimo.App()


@app.cell
def _():
    import marimo as mo
    return


@app.cell
def _():
    import pandas as pd
    import pybedtools as pbt
    return pbt, pd


@app.cell
def _(pbt):
    # Load gene tss annotations
    tss_bed = pbt.BedTool("../../annotations/gencode_tss.bed")
    genome_file = "../../../../Annotations/HG38/GRCh38_EBV.chrom.sizes.no.alt.tsv"

    regulatory_region_pybed = tss_bed.slop(
            g=genome_file,
            s=True,
            l=0,
            r=500
        ).sort(g=genome_file)

    # For each gene select the first promoter region only
    regulatory_region_df = regulatory_region_pybed.to_dataframe(names=[
        "chrom",
        "start",
        "end",  
        "transcript_id",
        "score",
        "strand",
        "gene_id",
        "gene_type",
        "gene_symbol"
    ])
    regulatory_region_df = regulatory_region_df.groupby("gene_id").first().reset_index()
    regulatory_region_df.head()
    return (regulatory_region_df,)


@app.cell
def _(pd):
    select_day = 4
    # Load discovery analysis results
    results_df = pd.read_csv(
        f"../../results/04-2025-10-27_sceptre_analysis/2025-11-13_sceptre_outputs/day{select_day}/sceptre_outputs/sceptre_discovery_results.csv",
        sep=",",
        header=0,
    )

    results_df["effect_size"] = (results_df["fold_change"] - 1) * 100
    # Keep only significant results
    results_df = results_df[results_df["significant"] == True]
    results_df.head()
    return results_df, select_day


@app.cell
def _(results_df):
    results_df.significant.value_counts()
    return


@app.cell
def _(regulatory_region_df, results_df, select_day):
    # BEDPE file format columns:
    # chr	start	end	gene_chr	gene_tss_start	gene_tss_500bp	eg_name	effect_size	day	color
    results_df["element_chrom"] = results_df["grna_target"].str.split(":").str[0]
    results_df["element_start"] = results_df["grna_target"].str.split(":").str[1].str.split("-").str[0].astype(int)-1
    results_df["element_end"] = results_df["grna_target"].str.split(":").str[1].str.split("-").str[1].astype(int)

    annotated_results = results_df.merge(
        regulatory_region_df[["chrom", "start", "end", "gene_id", "gene_symbol", "transcript_id"]],
        left_on="response_id",
        right_on="gene_symbol",
        how="left",
        suffixes=("", "_y")
    )

    annotated_results["day"] = 0
    annotated_results["pair_id"] = (
        annotated_results["gene_symbol"].astype(str)
        + "|" +
        annotated_results["grna_target"].astype(str)
    )

    annotated_results[[
        "element_chrom",
        "element_start",
        "element_end",
        "chrom",
        "start",
        "end",
        "pair_id",
        "gene_symbol",
        "transcript_id",
        "gene_id",
        "effect_size",
        "day",
        "significant"
    ]].to_csv(
        f"../../results/04-2025-10-27_sceptre_analysis/sceptre_discovery_results_day{select_day}_bedpe.tsv",
        sep="\t",
        header=True,
        index=False
    )
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
