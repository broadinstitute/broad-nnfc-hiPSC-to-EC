# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "marimo>=0.23.2",
#     "numpy>=2.4.4",
#     "pandas>=3.0.2",
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
    # Annotate SCEPTRE discovery results

    This notebook takes the output of SCEPTRE discovery and annotates it:
    - effect_size
    - se_effect_size
    - distance from targeted TSS
    - element genomic annotation
    - if inside a gene body, the gene name and distance to TSS
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
        element_annotations: Path = field(
            alias="--element-annotations",
            metadata={"help": "Path to element annotation file"}
        )
        quadrant_annotations: Path = field(
            alias="--quadrant-annotations",
            metadata={"help": "Path to quadrant annotation file"}
        )
        gene_trends_file: Path = field(
            alias="--gene-trends-file",
            metadata={"help": "Path to gene trends file"}
        )
        output_file: Path = field(
            default=Path("annotated_results.tsv"),
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
        --element-annotations element_annotations.tsv \\
        --gene-trends-file gene_trends.tsv \\
        --quadrant-annotations quadrant_annotations.tsv \\
        --output-file annotated_results.tsv


      # Use both TSS and element midpoints
      sceptre_results_to_bedpe \\
        --input-file discovery.tsv \\
        --tss-file tss_annotations.tsv \\
        --element-annotations element_annotations.tsv \\
        --quadrant-annotations quadrant_annotations.tsv \\ 
        --gene-trends-file gene_trends.tsv \\
        --output-file annotated_results.tsv \\
        --tss-midpoint \\
        --element-midpoint

    Notes:
    - Input file should be a SCEPTRE discovery results CSV file.
    - TSS file should contain genomic coordinates for transcription start sites.
    - Element annotations file should contain the following columns:
        - chrom
        - start 
        - end
        - element_id
        - annotation_type
        - gene_annotation
        - abs_distance_to_tss
    - Quadrant annotations file should contain the following columns:
        - ensembl_id
        - mean_upm
        - mean_umi
        - quadrant

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
                element_annotations=Path("../../results/03-2025-10-24_create_sceptre_input/elements.annotated.provenance.with_tss_distance.tsv"),
                quadrant_annotations=Path("../../results/2026-04-21_load_counts_and_qc/day0_ribo20_tumi55_umi1500/gene_umi_comparison_quadrants.tsv"),
                gene_trends_file=Path("../../results/gene_classification_by_Noam/gene_grouping_final.tsv"),
                output_file=Path("../../results/2026-04-22_sceptre_results/day4_grna20/sceptre_discovery_results_day4_grna20_annotated.tsv"), 
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
    import numpy as np
    import pandas as pd

    return np, pd


@app.cell
def _(args, pd):
    # Load gene TSS annotations
    tss_df = pd.read_csv(
        args.tss_file,
        sep="\t",
        header=0,
        names=[
            "tss_chrom",
            "tss_start",
            "tss_end",
            "gene_symbol",
            "score",
            "strand",
            "gene_id",
            "gene_type",
        ],
        dtype={
            "tss_start": "Int64",
            "tss_end": "Int64",
        },
    )


    # If tss_midpoint is True, use the midpoint of the TSS annotation as the TSS coordinate
    if args.tss_midpoint:
        mid = tss_df["tss_start"] + (tss_df["tss_end"] - tss_df["tss_start"]) // 2
        tss_df["tss_start"] = mid
        tss_df["tss_end"] = mid + 1

    # Select the first TSS for each gene (in case there are multiple TSS annotations per gene)
    tss_df = tss_df.groupby("gene_id", as_index=False).first()

    tss_df.head()
    return (tss_df,)


@app.cell
def _(args, pd):
    # Read elements annotations file
    element_annotations_df = pd.read_csv(args.element_annotations, sep="\t")
    element_annotations_df.head()
    # Create a dictionary mapping element_id to annotation info
    element_annotations_df = element_annotations_df.drop_duplicates(subset=["element_id"], keep="first")
    element_annotations_dict = element_annotations_df.set_index("element_id").to_dict(orient="index")
    # first few rows of the dictionary
    {element_id: {"chrom": row["chrom"], "start": row["start"], "end": row["end"], "annotation_type": row["annotation_type"], "gene_annotation": row["gene_annotation"], "abs_distance_to_tss": row["abs_distance_to_tss"]} for element_id, row in list(element_annotations_dict.items())[:3]}

    return (element_annotations_df,)


@app.cell
def _(args, pd):
    # Load quadrant annotations file
    quadrant_annotations_df = pd.read_csv(args.quadrant_annotations, sep="\t")
    quadrant_annotations_df.rename(columns={"ensembl_id": "gene_id", "mean_upm": "reference_upm", "mean_umi.x": "tap_umi"}, inplace=True)
    quadrant_annotations_df.head()
    return (quadrant_annotations_df,)


@app.cell
def _(args, pd):
    # Load gene trends file
    gene_trends_df = pd.read_csv(args.gene_trends_file, sep="\t")
    gene_trends_df.head()
    return (gene_trends_df,)


@app.cell
def _(
    args,
    element_annotations_df,
    gene_trends_df,
    np,
    pd,
    quadrant_annotations_df,
    tss_df,
):
    # Load SCEPTRE discovery results
    sceptre_results_df = pd.read_csv(args.input_file)

    # Compute effect size and standard error of effect size for each result
    sceptre_results_df["effect_size"] = (sceptre_results_df["fold_change"] - 1) * 100
    sceptre_results_df["se_effect_size"] = (sceptre_results_df["se_fold_change"] * 100)
    # Compute z-score for effect size
    sceptre_results_df["z_effect_size"] = sceptre_results_df["effect_size"] / sceptre_results_df["se_effect_size"]
    # Create a unique identifier for each target-gene pair
    sceptre_results_df["pair_id"] = sceptre_results_df["grna_target"] + "||" + sceptre_results_df["response_id"]

    # Merge with element annotations
    sceptre_results_df = sceptre_results_df.merge(
        element_annotations_df[["element_id", "chrom", "start", "end", "annotation_type", "gene_annotation", "abs_distance_to_tss"]],
        left_on="grna_target",
        right_on="element_id",
        how="left"
    )
    sceptre_results_df.rename(columns={"chrom": "element_chrom", "start": "element_start", "end": "element_end", "gene_annotation": "in_gene_annotation", "abs_distance_to_tss": "abs_distance_to_in_gene_tss"}, inplace=True)
    sceptre_results_df[["in_gene_id", "in_gene_symbol"]] = sceptre_results_df["in_gene_annotation"].str.split("|", expand=True)

    # Annotate link type. If response_id is identical to in_gene_symbol, self pair
    sceptre_results_df["is_self_pair"] = np.where(sceptre_results_df["response_id"] == sceptre_results_df["in_gene_symbol"], True, False)
    # Annotate whether the element is proximal to the TSS of the response gene (<= 1kb)
    sceptre_results_df["is_element_proximal"] = np.where(sceptre_results_df["abs_distance_to_in_gene_tss"] <= 1000, True, False)

    # Merge with TSS annotations to get TSS coordinates for each response gene
    sceptre_results_df = sceptre_results_df.merge(
        tss_df[["tss_chrom", "tss_start", "tss_end", "gene_id", "gene_symbol"]],
        left_on="response_id",
        right_on="gene_symbol",
        how="left",
        suffixes=("", "_y")
    )

    tss_pos = sceptre_results_df["tss_start"]

    if args.element_midpoint:
        elem_mid = (
            sceptre_results_df["element_start"]
            + (sceptre_results_df["element_end"] - sceptre_results_df["element_start"]) // 2
        )
        sceptre_results_df["e_g_abs_distance"] = (elem_mid - tss_pos).abs()
    else:
        sceptre_results_df["e_g_abs_distance"] = np.maximum.reduce([
            sceptre_results_df["element_start"] - tss_pos,
            tss_pos - sceptre_results_df["element_end"],
            0
        ])

    #Drop intermediate columns
    sceptre_results_df = sceptre_results_df.drop(columns=["tss_chrom", "tss_start", "tss_end", "gene_symbol"])
    # Rename columns for clarity
    sceptre_results_df.rename(columns={"response_id":"response_symbol", "gene_id": "response_id"}, inplace=True)

    # Merge with quadrant annotation
    sceptre_results_df = sceptre_results_df.merge(
        quadrant_annotations_df[["gene_id", "reference_upm", "tap_umi", "quadrant"]],
        left_on="response_id",
        right_on="gene_id",
        how="left"
    )

    # Merge with gene_trend
    sceptre_results_df = sceptre_results_df.merge(
        gene_trends_df[["gene", "group"]],
        left_on="response_symbol",
        right_on="gene",
        how="left"
    )

    # More column cleanup
    sceptre_results_df = sceptre_results_df.drop(columns=["gene_id", "gene"])
    # Rename columns for clarity
    sceptre_results_df.rename(columns={"group":"gene_trend"}, inplace=True)

    # Save to output file
    sceptre_results_df.to_csv(args.output_file, sep="\t", index=False)
    return


if __name__ == "__main__":
    app.run()
