# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "marimo>=0.19.10",
#     "pandas>=3.0.0",
#     "pybedtools>=0.12.0",
#     "pyzmq>=27.1.0",
# ]
# ///

import marimo

__generated_with = "0.19.11"
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
def _(pd):
    def genomic_distance(
        chrom1, start1, end1, chrom2, start2, end2, strand1=None, strand2=None
        ):
        """Return distance (bp) between two genomic intervals.

        - If intervals overlap, returns 0.
        - If on different chromosomes or missing coordinates, returns pd.NA.
        """
        if pd.isna(chrom1) or pd.isna(chrom2):
            return pd.NA
        if chrom1 != chrom2:
            return pd.NA
        if pd.isna(start1) or pd.isna(end1) or pd.isna(start2) or pd.isna(end2):
            return pd.NA
        start1 = int(start1)
        end1 = int(end1)
        start2 = int(start2)
        end2 = int(end2)
        if end1 < start2:
            return start2 - end1
        if end2 < start1:
            return start1 - end2
        return 0

    return (genomic_distance,)


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
    #results_fnp = "../../results/2026-01-22_create_sceptre_object_and_qc/background_permutations/day0_grna20/sceptre_discovery_results.csv"
    results_fnp = f"../../results/2026-01-22_create_sceptre_object_and_qc/day{select_day}_grna20/sceptre_discovery_results.csv"
    # Load discovery analysis results
    results_df = pd.read_csv(
        results_fnp,
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
def _(genomic_distance, regulatory_region_df, results_df, select_day):
    # BEDPE file format columns:
    # element_chrom	element_start	element_end	gene_chrom	gene_start	gene_end	pair_id	score	strand1	strand2	gene_symbol	transcript_id	gene_id	effect_size	day	    significant distance
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

    annotated_results["day"] = select_day
    annotated_results["pair_id"] = (
        annotated_results["gene_symbol"].astype(str)
        + "|" +
        annotated_results["grna_target"].astype(str)
    )

    # Use the absolute value of the effect size as the score for the BEDPE file (since BEDPE score is typically used for visualization and we want to visualize both positive and negative effects)
    annotated_results["score"] = annotated_results["effect_size"].abs()
    # Normalize betwen 0 and 1 using the max absolute effect size per day
    max_effect_size = annotated_results["effect_size"].abs().max()
    annotated_results["score"] = annotated_results["score"] / max_effect_size
    annotated_results["strand1"] = "."
    annotated_results["strand2"] = "."

    #distance from element to gene tss
    annotated_results["distance"] = annotated_results.apply(
        lambda row: genomic_distance(
            chrom1=row["element_chrom"],
            start1=row["element_start"],
            end1=row["element_end"],
            chrom2=row["chrom"],
            start2=row["start"],
            end2=row["end"],
            strand1=row["strand1"],
            strand2=row["strand2"]
        ),
        axis=1
    )


    annotated_results[[
        "element_chrom",
        "element_start",
        "element_end",
        "chrom",
        "start",
        "end",
        "pair_id",
        "score",
        "strand1",
        "strand2",
        "gene_symbol",
        "transcript_id",
        "gene_id",
        "effect_size",
        "day",
        "significant",
        "distance"
    ]].to_csv(
        f"../../results/2026-02-12_post_sceptre_analysis/sceptre_discovery_results_day{select_day}.bedpe",
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
