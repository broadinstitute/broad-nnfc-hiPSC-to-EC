# Using the 300bp input for design peaks
#
# These commands assume you're running from:
#   analyses/2026-04-14_refine_element_boundaries/
# and using bash (process substitution is used below).

# --- Parameters---
ROOT=../..
OUTDIR="$ROOT/results/2026-04-14_refine_element_boundaries"
GUIDES_BED="$ROOT/results/02-2025-10-21_align_guides/guides_uniquely_mapped.bed"
LOCI_BED="$ROOT/data/loci/ec_screen_16_loci_2mb_locus_final.bed"
DESIGN_PEAKS_DIR="$ROOT/data/design_peaks"
OLD_ELEMENT_ANN="$ROOT/results/03-2025-10-24_create_sceptre_input/all_element_annotations_final.bed"
MAX_GUIDE_DIFF=5

mkdir -p "$OUTDIR"

# Find all the guides falling in loci
bedtools intersect \
  -a "$GUIDES_BED" \
  -b "$LOCI_BED" \
  -wa \
> "$OUTDIR/guides_in_loci.bed"

# Merge all peaks and remove exact start and end
cat "$DESIGN_PEAKS_DIR"/* \
  | sort -k1,1 -k2,2n -u \
> "$OUTDIR/union_all_days_remove_exact_start_and_end.bed"

# Intersect the guides with the peaks
bedtools intersect \
  -a "$OUTDIR/union_all_days_remove_exact_start_and_end.bed" \
  -b "$GUIDES_BED" \
  -wa -wb \
> "$OUTDIR/peaks_x_guides.bed"

# Perform the smart union
# Merges overlapping peaks if the symmetric difference between guide sets is < MAX_GUIDE_DIFF
python merge_by_guides.py \
  "$OUTDIR/peaks_x_guides.bed" \
  "$OUTDIR/merged_by_guides.bed" 

# Guides that are inside loci but now missing from merged_by_guides
# This is likely caused by extra regions added for validation
bedtools intersect \
  -a "$OUTDIR/guides_in_loci.bed" \
  -b "$OUTDIR/merged_by_guides.bed" \
  -v \
> "$OUTDIR/missing_guides.bed"

# Taking the elements for these guides from the old annotations
bedtools intersect \
  -a "$OUTDIR/missing_guides.bed" \
  -b "$OLD_ELEMENT_ANN" \
  -wb \
  | cut -f9-13 \
  | sort -u \
> "$OUTDIR/missing_elements.bed"

# Get the guides overlapping the missing elements
bedtools intersect \
  -a <( cut -f1-3 "$OUTDIR/missing_elements.bed" ) \
  -b "$GUIDES_BED" \
  -wa -wb \
> "$OUTDIR/missing_elements_x_guides.bed"

cat "$OUTDIR/missing_elements_x_guides.bed" "$OUTDIR/peaks_x_guides.bed" \
  | sort -k1,1 -k2,2n \
> "$OUTDIR/merged_by_guides_with_missing_elements.bed"

# Final set for the analysis
python merge_by_guides.py \
  "$OUTDIR/merged_by_guides_with_missing_elements.bed" \
  "$OUTDIR/elements_bed_final.bed"
