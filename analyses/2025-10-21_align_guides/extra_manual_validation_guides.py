import marimo

__generated_with = "0.17.8"
app = marimo.App()


@app.cell
def _():
    import marimo as mo
    return


@app.cell
def _():
    import pandas as pd
    return (pd,)


@app.cell
def _(pd):
    # Original old guide metadata
    df = pd.read_csv('../../metadata/legacy/240918_EC_TAPseq_CRISPR.full.tsv',sep="\t")
    df.head()
    return (df,)


@app.cell
def _(pd):

    # Unique alignments
    mapping = pd.read_csv('../../results/02-2025-10-21_align_guides/ec_differentiation_with_pam_and_leading_g_hg38_filtered_all.bed', sep="\t", header=None, names=['chr','start','end', 'name', 'score','strand', 'NM', 'AS'])
    mapping.head()
    return (mapping,)


@app.cell
def _(pd):
    # Discarded alignments

    discarded_mapping = pd.read_csv('../../results/02-2025-10-21_align_guides/discarded.bed',sep="\t",header=None,names=['chr','start','end', 'name', 'score','strand'])
    discarded_mapping.name.nunique()
    return (discarded_mapping,)


@app.cell
def _(pd):
    # Unmapped guides
    unmapped_guides = pd.read_csv('../../results/02-2025-10-21_align_guides/unmapped.tsv', sep="\t", header=0)
    unmapped_guides.head()
    return (unmapped_guides,)


@app.cell
def _(discarded_mapping, mapping):
    # Find guide without a valid alignment
    discarded = set(discarded_mapping.name) - set(mapping.name)
    len(discarded)
    return (discarded,)


@app.cell
def _(df, unmapped_guides):
    # What was the annotation of the unmapped guides?
    merged_unmapped = unmapped_guides.merge(df, left_on='read_name', right_on='OligoID', how='left')
    merged_unmapped = merged_unmapped.drop_duplicates(subset=['read_name'])
    # Double check the ones that are not labeled as negative controls
    unexpected = merged_unmapped[~merged_unmapped['guideSet'].isin(['negative_control'])]
    unexpected
    return merged_unmapped, unexpected


@app.cell
def _(merged_unmapped):
    merged_unmapped[merged_unmapped['guideSet'].isin(['negative_control'])]
    return


@app.cell
def _(pd):
    # Check the blat results for the unmapped guides that are not a negative controls and the discarded ones
    psl_columns = [
        "matches","misMatches","repMatches","nCount",
        "qNumInsert","qBaseInsert","tNumInsert","tBaseInsert",
        "strand","qName","qSize","qStart","qEnd","tName",
        "tSize","tStart","tEnd","blockCount","blockSizes",
        "qStarts","tStarts"
    ]

    blat_results = pd.read_csv(
        '../../results/02-2025-10-21_align_guides/blat_complete_output.psl',
        sep=r"\s+",
        names=psl_columns,
        comment="#",
        skiprows=5,
        engine="python"
    )

    blat_results.head()
    return (blat_results,)


@app.cell
def _(discarded, unexpected):
    # Extract blat results for the unexpected unmapped guides and the discarded ones
    guides_to_check = set(unexpected['read_name']).union(discarded)
    len(guides_to_check)
    return


@app.cell
def _(blat_results, unexpected):
    # Double checking the unmapped
    blat_unmapped = blat_results[blat_results['qName'].isin(set(unexpected.read_name))]
    # Filter out all the alignment with insertions
    blat_unmapped = blat_unmapped[blat_unmapped['blockCount'] == 1]
    # Filter out all the alignments that do not cover the whole guide
    blat_unmapped = blat_unmapped[(blat_unmapped['qStart'] == 0) & (blat_unmapped['qEnd'] == blat_unmapped['qSize'])]
    # Filter out all the alignments that have any insertion or deletion
    blat_unmapped = blat_unmapped[(blat_unmapped['misMatches'] == 1) & (blat_unmapped['nCount'] == 1)]
    # This single one should be rescued
    blat_unmapped

    # BLAT results for the rescued unmapped guides:
    #chr17	38580078	38580101    240911_EC_TAPseq_Enhancers_31616
    # It does overlap a SINE
    # name: MIRc
    # location: chr2:23593987-23594095
    # score: 1000.0
    # Name: MIRc
    # Class: SINE
    # Family: MIR
    return


@app.cell
def _(blat_results, discarded):
    # Double checking the discarded
    blat_discarded = blat_results[blat_results['qName'].isin(set(discarded))]
    # Filter out all the alignment with insertions
    blat_discarded = blat_discarded[blat_discarded['blockCount'] == 1]
    # Filter out all the alignments that do not cover the whole guide
    blat_discarded = blat_discarded[(blat_discarded['qStart'] == 0) & (blat_discarded['qEnd'] == blat_discarded['qSize'])]
    # Filter out all the alignments that have any insertion or deletion
    blat_discarded = blat_discarded[(blat_discarded['misMatches'] == 1) & (blat_discarded['nCount'] == 1)]
    # Two can be rescued
    blat_discarded

    # BLAT results for the rescued discarded guides:
    # chr18	58104433	58104457    240911_EC_TAPseq_Enhancers_23658
    # chr17	38465058	38465081    240911_EC_TAPseq_Enhancers_31157
    # GEM mapper alignemts for the rescued discarded guide:
    # 240911_EC_TAPseq_Enhancers_31157	chr17	38465059	0	0	+	19M2I2M	2	13	21	discarded_protospacer_indel
    # 240911_EC_TAPseq_Enhancers_23658	chr18	58104434	256	0	+	21M2I1M	2	14	22	discarded_protospacer_indel
    return


@app.cell
def _(mapping):
    mapping.name.nunique()
    return


@app.cell
def _(blat_mapping_results, mapping):
    # BLAT does not report any alignment for 17 guides:
    set(mapping.name) - set(blat_mapping_results.qName)
    return


@app.cell
def _():
    # One of the example guide not mapper by BLAT:
    #chr7	44104369	44104389	240913_EC_TAPseq_TSS_849	0	+	1	18
    # BLAST confirms the location found by GEM
    return


@app.cell
def _(discarded, mapping, negative_control_guides):
    # 11 Negative controls are discarded. 451 + 11 = 462. That means that one negative control guide is mapped uniquely. 
    len(discarded.intersection(set(negative_control_guides['OligoID'])) - set(mapping.name) )
    return


@app.cell
def _(df, mapping):
    # At this point it is clear to me that BLAT is not sensitive enough to map all the guides, so I will rely on GEM results only
    # Only last check are the guides labeled as negative controls that map uniquely with GEM
    negative_control_guides = df[df['guideSet'] == 'negative_control']
    negative_control_mapped = mapping[mapping['name'].isin(set(negative_control_guides['OligoID']))]
    negative_control_mapped
    return (negative_control_guides,)


@app.cell
def _():
    # 240913_EC_TAPseq_TSS_5174 has an insertion where the N is.
    # Maybe we can remove it
    return


@app.cell
def _(mapping):
    # Keep only the guides with one valid alignment 43284 in total
    # 43158 with unique alignment.
    mapping_unique = mapping.groupby('name').filter(lambda x: len(x) == 1)
    mapping_unique.name.nunique()
    return (mapping_unique,)


@app.cell
def _(mapping_unique):
    # Checking the CIGAR strings of the unique mappings
    mapping_unique.NM.value_counts()

    return


@app.cell
def _():
    # In conclusion:
    # Rely only on GEM mapper results
    # Use the guides with no alignment as negative controls
    # Discard guides with multiple alignments
    # Keep only guides with one unique alignment with up to 2 mismatch as discovery pairs
    return


@app.cell
def _():
    # Deleted from HG19 to HG38:
    # 240911_EC_TAPseq_Enhancers_31605
    # 240911_EC_TAPseq_Enhancers_31606
    # 240911_EC_TAPseq_Enhancers_31607
    return


if __name__ == "__main__":
    app.run()
