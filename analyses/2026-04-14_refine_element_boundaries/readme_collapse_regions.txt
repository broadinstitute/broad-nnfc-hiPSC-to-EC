# Using the 300bp input for design peaks

# Find all the guides falling in loci
bedtools intersect -a ../02-2025-10-21_align_guides/guides_uniquely_mapped.bed -b ../../data/loci/ec_screen_16_loci_2mb_locus_final.bed -wa > guides_in_loci.bed

# Merge all peaks and remove exact start and end
cat ../../data/design_peaks/* | sort -k1,1 -k2,2n -u > union_all_days_remove_exact_start_and_end.bed

# Intersect the guides with the peaks
bedtools intersect \
  -a union_all_days_remove_exact_start_and_end.bed \
  -b ../02-2025-10-21_align_guides/guides_uniquely_mapped.bed \
  -wa -wb \
> peaks_x_guides.tsv

# Perform the smart union
# If guides(A) - guides(B) < 5 merge the peaks
python ../../analyses/2026-04-14_refine_element_boundaries/merge_by_guides.py peaks_x_guides.tsv merged_by_guides.bed

# Guides that are inside loci but now in the design peaks_x_guides
# This is likely caused by extra regions added for validation
bedtools intersect -a guides_in_loci.bed -b merged_by_guides.bed -v > missing_guides.bed 

# Taking the elements for these guides from the old annotations
bedtools intersect -a missing_guides.bed -b ../03-2025-10-24_create_sceptre_input/all_element_annotations_final.bed -wb | cut -f9-13 | sort -u > missing_elements.bed
# Final set for the analysis
cat missing_elements.bed <(cut -f1-5 merged_by_guides.bed) | sort -k1,1 -k2,2n | bedtools merge -i -  | awk -v OFS="\t" '{print $1,$2,$3,$1":"$2"-"$3}'> elements_bed_final.bed