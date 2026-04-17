OUTDIR=../../results/03-2025-10-24_create_sceptre_input

A=${OUTDIR}/elements.annotated.tsv
B=${OUTDIR}/all_element_annotations_final.bed

bedtools intersect -a "${A}" -b "${B}" -wa -wb | cut -f 1-6,11 | sort -k1,1 -k2,2n | sort -u > "${OUTDIR}/elements.annotated.provenance.tsv"