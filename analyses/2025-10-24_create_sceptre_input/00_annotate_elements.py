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
def _(np, pd):
    def sort_genomic(df, col='locus', place_non_genomic='last'):
        """
        Sort a DataFrame by genomic coordinates encoded as 'chrN:start-end'.
        Non-genomic labels (e.g., 'non-targeting') are placed at the end (default).
    
        Args:
            df: pandas DataFrame
            col: column name containing genomic coordinates as strings
            place_non_genomic: 'last' or 'first'
        
        Returns:
            Sorted DataFrame
        """

        # Extract only lines that match chr:start-end
        tmp = df[col].str.extract(r'^(chr[\w]+):(\d+)-(\d+)$')
        tmp.columns = ['chrom', 'start', 'end']

        # Convert numeric fields, NaN for non-genomic entries
        tmp['start'] = pd.to_numeric(tmp['start'], errors='coerce')
        tmp['end']   = pd.to_numeric(tmp['end'], errors='coerce')

        # Chromosome numeric sorting key
        def chrom_key(c):
            if pd.isna(c):
                return np.inf   # non-genomic → sorted last by default
            c = c.replace("chr", "")
            if c.isdigit():
                return int(c)
            return {'X': 23, 'Y': 24, 'M': 25, 'MT': 25}.get(c, 999)

        tmp['chrom_key'] = tmp['chrom'].apply(chrom_key)

        # OPTIONAL: place non-genomic entries first
        if place_non_genomic == 'first':
            tmp['chrom_key'] = tmp['chrom_key'].fillna(-1)
        else:
            tmp['chrom_key'] = tmp['chrom_key'].fillna(np.inf)

        # Sort by chromosomes → start → end
        df_sorted = (
            df.assign(
                chrom_key=tmp['chrom_key'],
                start=tmp['start'],
                end=tmp['end'],
            )
            .sort_values(['chrom_key', 'start', 'end'], na_position='last')
            .drop(columns=['chrom_key', 'start', 'end'])
        )

        return df_sorted
    return (sort_genomic,)


@app.cell
def _():
    # Mapped guides
    mapped_unique_guide_fnp = "../../metadata/uniquely_mapped_guides_final.bed"
    # Previously validated pairs file
    validated_pairs_fnp = "../../metadata/legacy/validate_pairs_from_previous_studies_final.tsv"
    # Previously validate promoters-targeting elements
    promoter_regions_fnp = "../../metadata/legacy/extra_list_of_targeted_elements_hg38_merged_final.bed"
    # Lifted over loci file
    loci_hg38_fnp = "../../metadata/legacy/ec_screen_16_loci_2mb_locus_final.bed"
    # Lifted over elements file
    lifted_elements_fnp = "../../metadata/legacy/240918_columns_10_13_26_noTSS_hg38_final_merged.bed"

    # Non-targeting control guides file
    non_targeting_guides_fnp = "../../results/02-2025-10-21_align_guides/non_targeting_guides_final.tsv"

    # Chromosome size file
    genome_file = "../../../../Annotations/HG38/GRCh38_EBV.chrom.sizes.no.alt.tsv"
    # TSS annotation file from GENCODE v44
    tss_file = "../../annotations/gencode_tss.bed"
    # Guide metadata file
    cellranger_feature_reference_file = "../../metadata/hiPSC-EC_feature-reference.csv"
    # Genes metadata file
    metadata_genes_fnp = "../../metadata/metadata_genes_final.csv"
    return (
        lifted_elements_fnp,
        loci_hg38_fnp,
        mapped_unique_guide_fnp,
        non_targeting_guides_fnp,
        promoter_regions_fnp,
        validated_pairs_fnp,
    )


@app.cell
def _(mo):
    mo.md(r"""
    ## Creating the guide to element dataframe
    First I will identify all the guides that overlaps the 44 previously validated elements. 26 distal elements and 18 proximal.
    """)
    return


@app.cell
def _(mapped_unique_guide_fnp, pd, pybedtools):
    # Load the uniquely mapped guides
    unique_df = pd.read_csv(mapped_unique_guide_fnp, sep="\t", header=0)
    unique_bed = pybedtools.BedTool.from_dataframe(unique_df)
    unique_bed.to_dataframe().head()
    return unique_bed, unique_df


@app.cell
def _(unique_df):
    unique_df.shape
    return


@app.cell
def _(mo):
    mo.md(r"""
    ### Overlap validated E-G pairs
    First I will overlap the uniquely mapped guides with 26 validated E-G pairs previously identified and included in this experiment.
    """)
    return


@app.cell
def _(pd, pybedtools, unique_bed, validated_pairs_fnp):
    e2g_validated = pd.read_csv(validated_pairs_fnp, sep="\t")
    coords = (
        e2g_validated['HG38']
            .str.extract(r'^(chr[^:]+):(\d+)-(\d+)$')
            .rename(columns={0: 'chrom', 1: 'start', 2: 'end'})
    )
    coords = coords.astype({'start': 'int64', 'end': 'int64'})
    e2g_validated = pd.concat([e2g_validated, coords], axis=1)
    # convert to 0-based BED start
    e2g_validated["start"] = e2g_validated["start"] - 1

    # build a proper BedTool from a DataFrame with 4 BED columns
    e2g_validated_bed = pybedtools.BedTool.from_dataframe(
        e2g_validated[['chrom', 'start', 'end', 'HG38', 'Response_ID', 'Source']].dropna()
    )

    # intersect unique guides with validated elements
    e2g_validated_bedtools = (
        unique_bed
            .intersect(
                e2g_validated_bed,
                wa=True,
                wb=True
            )
    )

    e2g_validated_bedtools_df = e2g_validated_bedtools.to_dataframe(
        names=[
            "guide_chromosome",
            "guide_start",
            "guide_end",
            "guide_id",
            "score",
            "strand",
            "NM",
            "AS",
            "guide_alias",
            "element_chromosome",
            "element_start",
            "element_end",
            "element_name",
            "response_id",
            "source"
            ]
        )

    # number of unique guides that overlap any validated element
    overlapping_guide_ids = set(e2g_validated_bedtools_df['guide_id'])
    e2g_validated_bedtools_df
    return e2g_validated_bedtools_df, overlapping_guide_ids


@app.cell
def _(e2g_validated_bedtools_df, unique_df):
    f"{e2g_validated_bedtools_df.shape[0]} guides overlapping previously validated E-G pairs",f"{unique_df.shape[0]-e2g_validated_bedtools_df.shape[0]} remaining."
    return


@app.cell
def _(mo):
    mo.md(r"""
    ### Overlap with TSSs
    There are 18 promoters previously validated, but all the TSS can be considered positive controls.
    For this reason I am overlapping once against all the TSSs. The previously validated will be in the form:
    ```
    chr15:41221280-41221780_DLL4_TSS # previously validated
    THAP5_TSS                        # other TSSs
    ```

    WARNING: Guides have been designed to target TSSs of non-readout genes.
    If gene is in readout then it should be considered a positive control pair.
    """)
    return


@app.cell
def _(overlapping_guide_ids, pybedtools, unique_df):
    # Get the guides not overlapping validated E-G pairs (by guide_id, not by intervals)
    # overlapping_guide_ids was computed in the previous cell
    unique_df_step2 = unique_df[
        ~unique_df['guide_id'].isin(overlapping_guide_ids)
    ].reset_index(drop=True)

    unique_bed_step2 = pybedtools.BedTool.from_dataframe(unique_df_step2)

    unique_bed_step2.to_dataframe().shape
    return unique_bed_step2, unique_df_step2


@app.cell
def _(csv, pd, promoter_regions_fnp, pybedtools, unique_bed_step2):
    # Load previously validated promoters
    promoters_with_validated_guides = ['CDH5', 'CTNND1', 'DLL4', 'FAM136A', 'FJX1', 'HAND1', 'HEY2', 'KDR', 'MOGS', 'PRKD1', 'PYCARD', 'RBMXL1', 'SAT2', 'SOX2', 'SRPK1', 'TMEM107', 'TMEM256', 'VKORC1']
    promoter_regions_df = pd.read_csv(
        promoter_regions_fnp,
        sep="\t",
        names=[
            "chrom", "start", "end",
            "element_id", "element_list",
            "element_alias_list", "source"
        ],
        engine="python",
        quoting=csv.QUOTE_NONE,
    )

    promoter_regions_bed = pybedtools.BedTool.from_dataframe(promoter_regions_df)

    # intersect unique guides with validated elements
    promoter_validated_bedtools = (
        unique_bed_step2
            .intersect(
                promoter_regions_bed,
                wa=True,
                wb=True
            )
    )

    promoter_validated_bedtools_df = promoter_validated_bedtools.to_dataframe(
        names=[
            "guide_chromosome",
            "guide_start",
            "guide_end",
            "guide_id",
            "score",
            "strand",
            "NM",
            "AS",
            "guide_alias",
            "element_chromosome",
            "element_start",
            "element_end",
            "element_name",
            "element_list",
            "element_alias_list",
            "source"
            ]
        )

    promoter_validated_bedtools_df["is_validated"] = promoter_validated_bedtools_df['element_alias_list'].str.startswith('chr')
    promoter_validated_bedtools_df["response_id"] = promoter_validated_bedtools_df['element_alias_list'].str.replace("_TSS","")
    promoter_validated_bedtools_df.loc[promoter_validated_bedtools_df["is_validated"], "response_id"] = promoter_validated_bedtools_df['element_alias_list'][promoter_validated_bedtools_df["is_validated"]].str.split("_").str[1]
    promoter_validated_bedtools_df
    return (promoter_validated_bedtools_df,)


@app.cell
def _(promoter_validated_bedtools_df, unique_df_step2):
    f"{len(set(promoter_validated_bedtools_df.guide_id))} guides overlapping previously validated promoter pairs",f"{unique_df_step2.shape[0] - len(set(promoter_validated_bedtools_df.guide_id))} remaining."
    return


@app.cell
def _(promoter_validated_bedtools_df, pybedtools, unique_df_step2):
    unique_df_step3 = unique_df_step2[
        ~unique_df_step2['guide_id'].isin(set(promoter_validated_bedtools_df['guide_id']))
    ].reset_index(drop=True)

    unique_bed_step3 = pybedtools.BedTool.from_dataframe(unique_df_step3)

    unique_df_step3.shape
    return unique_bed_step3, unique_df_step3


@app.cell
def _(mo):
    mo.md(r"""
    ### Overlap with lifted over elements
    """)
    return


@app.cell
def _(lifted_elements_fnp, pybedtools, unique_bed_step3, unique_df_step3):
    lifted_elements_bed = pybedtools.BedTool(lifted_elements_fnp)

    lifted_elements_df = unique_bed_step3.intersect(lifted_elements_bed,wa=True,wb=True).to_dataframe(names=[
            "guide_chromosome",
            "guide_start",
            "guide_end",
            "guide_id",
            "score",
            "strand",
            "NM",
            "AS",
            "guide_alias",
            "element_chromosome",
            "element_start",
            "element_end",
            ])

    # build 'chr:(start+1)-end' element name (start from 0-based BED -> 1-based coord)
    lifted_elements_df["element_name"] = (
        lifted_elements_df["element_chromosome"].astype(str)
        + ":"
        + (lifted_elements_df["element_start"] + 1).astype(str)
        + "-"
        + lifted_elements_df["element_end"].astype(str)
    )

    unique_df_step4 = unique_df_step3[
        ~unique_df_step3['guide_id'].isin(set(lifted_elements_df['guide_id']))
    ].reset_index(drop=True)


    lifted_elements_df
    return lifted_elements_df, unique_df_step4


@app.cell
def _(lifted_elements_df, unique_df_step3):
    f"{len(set(lifted_elements_df.guide_id))} guides overlapping lifted over element regions",f"{unique_df_step3.shape[0] - len(set(lifted_elements_df.guide_id))} remaining."
    return


@app.cell
def _(loci_hg38_fnp, pybedtools, unique_df_step4):
    loci_bed = pybedtools.BedTool(loci_hg38_fnp)

    unique_bed_step4 = pybedtools.BedTool.from_dataframe(unique_df_step4)

    lifted_loci_elements_df = unique_bed_step4.intersect(loci_bed,wa=True,wb=True).to_dataframe(names=[
            "guide_chromosome",
            "guide_start",
            "guide_end",
            "guide_id",
            "score",
            "strand",
            "NM",
            "AS",
            "guide_alias",
            "element_chromosome",
            "element_start",
            "element_end",
            "element_name"
            ])


    # lifted_loci_elements_df should be empty so the guides in unique_df_step4 are the putative safe targeting guides
    assert lifted_loci_elements_df.shape[0] == 0
    unique_df_step4
    putative_safe_targeting_df = unique_df_step4.copy()
    putative_safe_targeting_df["element_name"] = "safe-targeting"
    putative_safe_targeting_df.rename(columns={"alias": "guide_alias"}, inplace=True)
    putative_safe_targeting_df
    return (putative_safe_targeting_df,)


@app.cell
def _(mo):
    mo.md(r"""
    ### Writing file to disk
    """)
    return


@app.cell
def _(
    e2g_validated_bedtools_df,
    lifted_elements_df,
    promoter_validated_bedtools_df,
    putative_safe_targeting_df,
):
    # Safe targeting guides
    putative_safe_targeting_df.to_csv("../../results/03-2025-10-24_create_sceptre_input/putative_safe_targeting_guides_final.tsv", sep="\t", index=False, header=True)
    # Lifted elements overlapping guides
    lifted_elements_df.to_csv("../../results/03-2025-10-24_create_sceptre_input/lifted_elements_overlapping_guides_final.tsv", sep="\t", index=False, header=True)
    # Promoters overlapping guides
    promoter_validated_bedtools_df.to_csv("../../results/03-2025-10-24_create_sceptre_input/promoters_overlapping_guides_final.tsv", sep="\t",  index=False, header=True)
    # E-G pairs overlapping guides
    e2g_validated_bedtools_df.to_csv("../../results/03-2025-10-24_create_sceptre_input/e2g_validated_overlapping_guides_final.tsv", sep="\t", index=False, header=True)
    return


@app.cell
def _(
    e2g_validated_bedtools_df,
    lifted_elements_df,
    non_targeting_guides_fnp,
    pd,
    promoter_validated_bedtools_df,
    putative_safe_targeting_df,
    sort_genomic,
):
    non_targeting_guides_df = pd.read_csv(non_targeting_guides_fnp, sep="\t", header=0)
    non_targeting_guides_df.rename(columns={"alias": "guide_alias"}, inplace=True)
    non_targeting_guides_df["element_name"] = "non-targeting"
    # Create on file with all the guide_alias and their corresponding element_name
    all_guides_annotations_df = pd.concat([
        putative_safe_targeting_df[["guide_id", "guide_alias", "element_name"]].assign(source="safe-targeting"),
        lifted_elements_df[["guide_id", "guide_alias", "element_name"]].assign(source="lifted_element"),
        promoter_validated_bedtools_df[["guide_id", "guide_alias", "element_name"]].assign(source="promoter"),
        e2g_validated_bedtools_df[["guide_id", "guide_alias", "element_name"]].assign(source="e2g_validated"),
        non_targeting_guides_df[["guide_id", "guide_alias", "element_name"]].assign(source="non-targeting")
        ])


    # Reorder columns and changing names because of SCEPTRE input requirements
    # SCEPTRE doesn't like 'safe-targeting' but I prefer not to remove it at this point
    all_guides_annotations_df = all_guides_annotations_df[["guide_alias", "element_name","guide_id","source"]]
    all_guides_annotations_df  = all_guides_annotations_df.rename(columns={"guide_alias": "grna_id", "element_name":"grna_target", "guide_id":"grna_sequence"})
    all_guides_annotations_df.sort_values(by=["grna_target"], inplace=True)
    all_guides_annotations_df = sort_genomic(all_guides_annotations_df, col="grna_target")
    all_guides_annotations_df.to_csv("../../results/03-2025-10-24_create_sceptre_input/all_guides_annotations_final.tsv", sep="\t", index=False, header=True)
    all_guides_annotations_df
    return


@app.cell
def _(
    e2g_validated_bedtools_df,
    lifted_elements_df,
    pd,
    promoter_validated_bedtools_df,
):
    # Writing the final BED file with the total gRNA_targets
    all_element_annotations_df = pd.concat([
        lifted_elements_df[["element_chromosome", "element_start", "element_end", "element_name"]].assign(source="lifted_element"),
        promoter_validated_bedtools_df[["element_chromosome", "element_start", "element_end", "element_name"]].assign(source="promoter"),
        e2g_validated_bedtools_df[["element_chromosome", "element_start", "element_end", "element_name"]].assign(source="e2g_validated")
        ]).drop_duplicates()
    all_element_annotations_df.to_csv("../../results/03-2025-10-24_create_sceptre_input/all_element_annotations_final.bed", sep="\t", index=False, header=False)

    all_element_annotations_df
    return


@app.cell
def _(
    alt,
    e2g_validated_bedtools_df,
    lifted_elements_df,
    pd,
    promoter_validated_bedtools_df,
):
    # Visualizing the distribution og guides per element
    all_elements = pd.concat([
        lifted_elements_df[["element_name"]].assign(source="lifted_element"),
        promoter_validated_bedtools_df[["element_name"]].assign(source="promoter"),
        e2g_validated_bedtools_df[["element_name"]].assign(source="e2g_validated"),
    ], ignore_index=True)

    # Get counts per element & source
    counts = (
        all_elements
        .groupby(["element_name", "source"])
        .size()
        .reset_index(name="count")
    )

    chart = (
        alt.Chart(counts)
        .transform_density(
            "count",
            groupby=["source"],
            as_=["count", "density"]
        )
        .mark_area(opacity=0.4)
        .encode(
            x="count:Q",
            y="density:Q",
            color="source:N"
        ).facet(column="source:N")
    )

    chart
    return (counts,)


@app.cell
def _(alt, counts):
    cdf = (
        alt.Chart(counts)
        .transform_window(
            sort=[{"field": "count"}],
            groupby=["source"],
            cumulative_count="count()"
        )
        .transform_joinaggregate(
            total="count()",
            groupby=["source"]
        )
        .transform_calculate(
            cdf="datum.cumulative_count / datum.total"
        )
        .mark_line()
        .encode(
            x=alt.X("count:Q", title="Occurrences"),
            y=alt.Y("cdf:Q", title="CDF (fraction ≤ x)", scale=alt.Scale(domain=[0, 1])),
            color="source:N"
        )
    )

    cdf
    return


if __name__ == "__main__":
    app.run()
