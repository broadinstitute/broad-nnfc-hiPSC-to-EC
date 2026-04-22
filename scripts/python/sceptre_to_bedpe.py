# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "marimo>=0.23.2",
#     "pandas>=3.0.0",
#     "pybedtools>=0.12.0",
#     "pyzmq>=27.1.0",
#     "simple-parsing>=0.1.8",
# ]
# ///

import marimo

__generated_with = "0.23.2"
app = marimo.App()


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell
def _():
    import simple_parsing

    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Convert a TSV file with SCEPTRE discovery results into a BEDPE file for visualization in genome browsers.
    """)
    return


@app.cell
def _():
    from dataclasses import dataclass
    from pathlib import Path
    from simple_parsing import field

    @dataclass
    class Args:
        input_file: Path = field(
            alias="--input-file",
            metadata={"help": "Path to SCEPTRE discovery results CSV file"}
        )
        tss_file: Path = field(
            alias="--tss-file",
            metadata={"help": "Path to TSS annotation file"}
        )
        genome_file: Path = field(
            alias="--genome-file",
            metadata={"help": "Path to genome file with chromosome sizes"}
        )
        output_file: Path = field(
            default=Path("output.bedpe"),
            alias="--output-file",
            metadata={"help": "Output BEDPE file path"}
        )
        tss_midpoint: bool = field(
            default=False,
            alias="--tss-midpoint",
            metadata={"help": "Use midpoint of the TSS annotation when computing distances"}
        )
        element_midpoint: bool = field(
            default=False,
            alias="--element-midpoint",
            metadata={"help": "Use midpoint of the element annotation when computing distances to TSS"}
        )

    epilog = """
    Examples:

      # Basic usage
      sceptre_results_to_bedpe \\
        --input-file discovery.tsv \\
        --tss-file tss_annotations.tsv \\
        --genome-file genome_file.tsv

      # Use both TSS and element midpoints
      sceptre_results_to_bedpe \\
        --input-file discovery.tsv \\
        --tss-file tss_annotations.tsv \\
        --genome-file genome_file.tsv \\
        --tss-midpoint \\
        --element-midpoint

    Notes:
      - Input file should be a SCEPTRE discovery results CSV file.
      - TSS file should contain genomic coordinates for transcription start sites.
      - Output is written in BEDPE format for downstream genomic analysis.
    """
    return Args, Path, epilog


@app.cell
def _(Args, Path, epilog, mo):
    from simple_parsing import parse

    def parse_args():
        if mo.running_in_notebook():
            return Args(
                # Input files for the notebook
                input_file=Path("../../results/2026-04-22_sceptre_results/day4_grna20/sceptre_discovery_results.csv"), 
                tss_file=Path("../../annotations/genes/gencode.v43.protein_coding.TSS500bp.bed"),
                genome_file=Path("../../annotations/GRCh38_EBV.chrom.sizes.no.alt.tsv"),
                output_file=Path("../../results/2026-04-22_sceptre_results/day4_grna20/sceptre_discovery_results_day4_grna20.bedpe"), 
                tss_midpoint=True,
                element_midpoint=True,
            )
        else:
            return parse(
                        Args,
                        description="Convert SCEPTRE discovery results to BEDPE format.",
                        epilog=epilog,
                    )


    return (parse_args,)


@app.cell
def _(mo, parse_args):
    args = parse_args()
    mo.json(vars(args))
    return (args,)


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

    def to_midpoint(interval):
        # Convert a genomic interval to its midpoint (single base pair).
        mid = interval.start + (interval.end - interval.start) // 2
        interval.start = mid
        interval.end = mid + 1  # BED is half-open
        return interval

    return genomic_distance, to_midpoint


@app.cell
def _(args, pbt, to_midpoint):
    # Load gene tss annotations
    tss_pybed = pbt.BedTool(args.tss_file)

    # If tss_midpoint is True, use the midpoint of the TSS annotation as the TSS coordinate
    if args.tss_midpoint:
        tss_pybed = tss_pybed.each(to_midpoint).saveas()


    tss_df = tss_pybed.to_dataframe(names=[
        "chrom",
        "start",
        "end",  
        "gene_symbol",
        "score",
        "strand",
        "gene_id",
        "gene_type",
    ])

    # Select the first TSS for each gene(in case there are multiple TSS annotations per gene)
    tss_df = tss_df.groupby("gene_id").first().reset_index()
    tss_df.head()
    return (tss_df,)


@app.cell
def _(args, pd):
    # Load discovery analysis results
    results_df = pd.read_csv(
        args.input_file,
        header=0,
    )

    # Compute effect size and standard error of effect size for each result
    results_df["effect_size"] = (results_df["fold_change"] - 1) * 100
    results_df["se_effect_size"] = (results_df["se_fold_change"] * 100)
    # Compute z-score for effect size
    results_df["z_effect_size"] = results_df["effect_size"] / results_df["se_effect_size"]
    # Create a unique identifier for each target-gene pair
    results_df["pair_id"] = results_df["grna_target"] + "||" + results_df["response_id"]

    # Keep only significant results
    results_df = results_df[results_df["significant"] == True]
    results_df.head()
    return (results_df,)


@app.cell
def _(results_df):
    results_df.significant.value_counts()
    return


@app.cell
def _(args, genomic_distance, results_df, tss_df):
    # BEDPE file format columns:
    # element_chrom	element_start	element_end	gene_chrom	gene_start	gene_end	pair_id	score	strand1	strand2	gene_symbol	transcript_id	gene_id	effect_size	day	    significant distance
    results_df["element_chrom"] = results_df["grna_target"].str.split(":").str[0]
    results_df["element_start"] = results_df["grna_target"].str.split(":").str[1].str.split("-").str[0].astype(int)-1
    results_df["element_end"] = results_df["grna_target"].str.split(":").str[1].str.split("-").str[1].astype(int)

    annotated_results = results_df.merge(
        tss_df[["chrom", "start", "end", "gene_id", "gene_symbol"]],
        left_on="response_id",
        right_on="gene_symbol",
        how="left",
        suffixes=("", "_y")
    )

    annotated_results["pair_id"] = (
        annotated_results["gene_symbol"].astype(str)
        + "|" +
        annotated_results["grna_target"].astype(str)
    )

    # Use the absolute value of the effect size as the score for the BEDPE file (since BEDPE score is typically used for visualization and we want to visualize both positive and negative effects)
    # Normalize betwen 0 and 1 using the max absolute effect size per day
    max_effect_size = annotated_results["effect_size"].abs().max()
    annotated_results["score"] = annotated_results["effect_size"].abs() / max_effect_size
    annotated_results["strand1"] = "."
    annotated_results["strand2"] = "."

    # Distance from element to gene tss
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

    # Column order for BEDPE file
    cols = [
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
        "gene_id",
        "effect_size",
        "significant",
        "distance"
    ]

    annotated_results[cols].to_csv(
        args.output_file,
        sep="\t",
        index=False,
        header=["#element_chrom", *cols[1:]],
    )
    return


if __name__ == "__main__":
    app.run()
