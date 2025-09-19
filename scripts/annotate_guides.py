# This script takes in input the results from the guide RNA alignment,
# merges together all the guides withing a distance threshold,
# and annotates the merged guides with the gene they target.
# If a guide does not have an alignment on the genome is annotated as 
# 'non-targeting'. If it aligns on the genome and it is around 2MB of
# a gene it is annotated as 'targeting', 'safe-targeting' otherwise.

import pandas as pd
import re
import pybedtools
import pysam

# Parameters
SAM_FILE = "results/2025-09-15/all_guides_hg38.sam"      # Input guide alignment SAM file
#   - SAM_FILE: Path to the SAM file containing guide RNA alignments.
#     Each record should represent a guide RNA mapped to the genome.
GENE_ANNOT_FILE = "results/2025-09-15/gencode_tss.bed" # Gene annotation BED file
#   - GENE_ANNOT_FILE: Path to a BED file with gene annotations.
#     Expected columns: chr, start, end, ..., gene_symbol (column 9).
DIST_THRESHOLD = 2000                         # Distance threshold for merging (bp)
#   - DIST_THRESHOLD: Maximum distance (in bp) between guides to merge into a single interval.
TARGET_WINDOW = 2_000_000                     # Window for 'targeting' annotation (bp)
#   - TARGET_WINDOW: Maximum distance (in bp) from a merged guide to a gene to be considered 'targeting'.

def parse_sam(sam_path: str) -> pd.DataFrame:
    """
    Parse a SAM file and extract guide alignments.

    Args:
        sam_path (str): Path to SAM file.

    Returns:
        pd.DataFrame: DataFrame with columns ['guide_id', 'chr', 'start', 'end'].
    """
    rows = []
    with pysam.AlignmentFile(sam_path, "r") as samfile:
        for aln in samfile:
            guide_id = aln.query_name
            if aln.is_unmapped:
                rows.append({'guide_id': guide_id, 'chr': None, 'start': None, 'end': None})
                continue
            chr = samfile.get_reference_name(aln.reference_id)
            # pysam: reference_start is 0-based, reference_end is 0-based exclusive
            start = aln.reference_start      # BED format expects 0-based start
            end = aln.reference_end          # BED format expects 0-based exclusive end
            rows.append({'guide_id': guide_id, 'chr': chr, 'start': start, 'end': end})
    return pd.DataFrame(rows)

def load_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load guide alignments and gene annotations.

    Returns:
        tuple: (guides_df, genes_df)
    """
    guides = parse_sam(SAM_FILE)
    # BED: chr, start, end, ..., gene_symbol (column 9)
    genes = pd.read_csv(GENE_ANNOT_FILE, sep='\t', header=None,
                        names=['chr', 'start', 'end', 'transcript_id', 'score', 'strand', 'gene_id', 'gene_type', 'gene_symbol'])
    return guides, genes

def merge_guides(guides: pd.DataFrame, dist_threshold: int) -> pd.DataFrame:
    """
    Merge guides within a distance threshold using pybedtools.

    Args:
        guides (pd.DataFrame): Guide intervals.
        dist_threshold (int): Merge distance in bp.

    Returns:
        pd.DataFrame: Merged intervals with guide IDs.
    """
    # Prepare BED format: chr, start, end, guide_id
    bed_df = guides.dropna(subset=['chr', 'start', 'end']).copy()
    # Ensure start/end are integers
    bed_df['start'] = bed_df['start'].astype(int)
    bed_df['end'] = bed_df['end'].astype(int)
    bed_df = bed_df[['chr', 'start', 'end', 'guide_id']]
    bedtool = pybedtools.BedTool.from_dataframe(bed_df)
    merged = bedtool.sort().merge(d=dist_threshold, c=4, o='distinct')
    merged_df = merged.to_dataframe(names=['chr', 'start', 'end', 'guide_ids'])
    # Add back unmapped guides
    unmapped = guides[guides['chr'].isna()].copy()
    unmapped['guide_ids'] = unmapped['guide_id']
    merged_df = pd.concat([merged_df, unmapped[['chr', 'start', 'end', 'guide_ids']]], ignore_index=True)
    return merged_df

def annotate_guides(merged_guides: pd.DataFrame, genes: pd.DataFrame, target_window: int) -> pd.DataFrame:
    """
    Annotate merged guides with gene targets.

    Args:
        merged_guides (pd.DataFrame): Merged guide intervals.
        genes (pd.DataFrame): Gene annotation intervals.
        target_window (int): Window for 'targeting' annotation.

    Returns:
        pd.DataFrame: Annotated guides.
    """
    def find_annotation(row):
        # If no alignment, non-targeting
        if pd.isna(row['chr']):
            return 'non-targeting', None
        # Find closest gene
        gene_hits = genes[genes['chr'] == row['chr']]
        gene_hits['dist'] = gene_hits.apply(
            lambda g: min(abs(row['start'] - g['start']), abs(row['end'] - g['end'])), axis=1
        )
        closest = gene_hits.loc[gene_hits['dist'].idxmin()] if not gene_hits.empty else None
        if closest is not None and closest['dist'] <= target_window:
            return 'targeting', closest['gene_symbol']
        else:
            return 'safe-targeting', None

    merged_guides[['annotation', 'target_gene']] = merged_guides.apply(
        lambda row: pd.Series(find_annotation(row)), axis=1
    )
    return merged_guides


def main() -> None:
    """
    Main entry point for guide annotation workflow.
    """
    guides, genes = load_data()
    merged_guides = merge_guides(guides, DIST_THRESHOLD)
    annotated_guides = annotate_guides(merged_guides, genes, TARGET_WINDOW)
    annotated_guides.to_csv("results/2025-09-15/annotated_guides.csv", index=False)
    print("Annotated guides saved to results/2025-09-15/annotated_guides.csv")

if __name__ == "__main__":
    main()
