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
    alt.data_transformers.enable("vegafusion")
    return alt, pd


@app.cell
def _(pd):
    # Import the metadata
    metadata_guide = pd.read_csv("metadata/hiPSC-EC/metadata_guides_with_accession.csv")
    metadata_guide["is_safe_targeting"] = metadata_guide["guideSet"] == "safe_targeting"
    metadata_guide["is_negative_control"] = metadata_guide["guideSet"] == "negative_control"

    metadata_guide
    return (metadata_guide,)


@app.cell
def _(pd):
    mapped = pd.read_csv("results/2025-10-01_ec_guide_annotation/ec_differentiation_with_pam_and_leading_g_hg38.bed", sep="\t")
    mapped
    return (mapped,)


@app.cell
def _(pd):
    # Alignment results
    unmapped = pd.read_csv("results/2025-10-01_ec_guide_annotation/ec_differentiation_with_pam_and_leading_g_hg38_unmapped.tsv", sep="\t")
    unmapped
    return (unmapped,)


@app.cell
def _(metadata_guide, mo, pd):
    def one_flag_table(df, col):
        vc = df[col].value_counts().reindex([True, False], fill_value=0)
        out = pd.DataFrame({"Value": ["True","False"], "Count": [vc.get(True,0), vc.get(False,0)]})
        out["%"] = (out["Count"] / out["Count"].sum()).map(lambda x: f"{x:.1%}")
        return out

    safe_tbl = one_flag_table(metadata_guide, "is_safe_targeting")
    neg_tbl  = one_flag_table(metadata_guide, "is_negative_control")

    mo.hstack(
        [mo.vstack(
            [
             mo.md("### is_safe_targeting"),
             mo.ui.table(safe_tbl)
            ]
        ),
         mo.vstack(
             [
              mo.md("### is_negative_control"),
              mo.ui.table(neg_tbl)
             ]
         )
        ]
    )
    return


@app.cell
def _(metadata_guide, unmapped):
    # Find the annotation of the unmapped guides
    metadata_guide["is_unmapped"] = metadata_guide["name"].isin(unmapped["guide_id"])
    # Check how many unmapped guides are safe targeting or negative control

    metadata_guide["guideSet"][metadata_guide["is_unmapped"]].value_counts()
    return


@app.cell
def _(mapped, metadata_guide):
    negative_controls_mapping = metadata_guide["name"][~metadata_guide["is_unmapped"] & metadata_guide["is_negative_control"]]

    # Look where the negative controls map on the genome
    negative_controls_mapping_locations = mapped[mapped["guide_id"].isin(negative_controls_mapping)]
    negative_controls_mapping_locations
    return


@app.cell
def _(alt, mapped):
    # Plot distribution of NM scores
    chart = alt.Chart(mapped).mark_bar().encode(
        x=alt.X("NM", bin=alt.Bin(maxbins=30), title="Number of Mismatches (NM)"),
        y=alt.Y("count()", title="Count of Guides",scale=alt.Scale(type="log")),
        tooltip=[alt.Tooltip("count()", title="Count of Guides")]
    ).properties(
        width=400,
        height=300,
        title="Distribution of NM scores"
    )
    chart.save('plots/distribution_NM_scores.png', scale_factor=4)

    chart
    return


@app.cell
def _(mapped):
    # Find all the guides that don't have any alternate alignments
    unique_alignments = mapped.groupby("guide_id").size().reset_index(name="alignment_count").query("alignment_count == 1")

    mapped[mapped["guide_id"].isin(unique_alignments["guide_id"])]["NM"].value_counts()
    return


@app.cell
def _(mapped):
    multiple_alignments = mapped.groupby("guide_id").size().reset_index(name="alignment_count").query("alignment_count > 1")
    mapped[mapped["guide_id"].isin(multiple_alignments["guide_id"])]
    return


@app.cell
def _(mapped, pd, unmapped):
    # Keep only primary alignments with NM <= 2
    # Check the guide ids of the discarded and if no alignmnets pass the filter add them to the unmapped
    mask_primary = mapped["type"] == "primary"
    mask_nm2 = mapped["NM"] <= 2
    filtered_mapped = mapped[mask_primary & mask_nm2]
    filtered_mapped
    discarded_guides = set(mapped["guide_id"]) - set(filtered_mapped["guide_id"])
    discarded_guides
    filtered_mapped
    # Add the discarded guides to the unmapped table
    new_unmapped = pd.DataFrame({"guide_id": list(discarded_guides)})
    new_unmapped["reason"] = "no_primary_with_NM<=2"
    new_unmapped
    final_unmapped = pd.concat([unmapped, new_unmapped], ignore_index=True)
    final_unmapped
    final_unmapped["guide_id"].nunique()
    # Save the final unmapped table
    final_unmapped.to_csv("results/2025-10-01_ec_guide_annotation/ec_differentiation_with_pam_and_leading_g_hg38_final_unmapped.tsv", sep="\t", index=False)
    # Save the filtered mapped table
    filtered_mapped.to_csv("results/2025-10-01_ec_guide_annotation/ec_differentiation_with_pam_and_leading_g_hg38_filtered.bed", sep="\t", index=False, header=False)
    return filtered_mapped, final_unmapped


@app.cell
def _(filtered_mapped, final_unmapped):
    filtered_mapped.shape,final_unmapped.shape
    return


@app.cell
def _(final_unmapped, metadata_guide):
    # Find the annotation of the unmapped guides
    metadata_guide["is_unmapped_extended"] = metadata_guide["name"].isin(final_unmapped["guide_id"])
    # Check how many unmapped guides are safe targeting or negative control

    metadata_guide["guideSet"][metadata_guide["is_unmapped_extended"]].value_counts()
    return


if __name__ == "__main__":
    app.run()
