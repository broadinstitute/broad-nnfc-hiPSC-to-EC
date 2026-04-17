OUTDIR=../../results/03-2025-10-24_create_sceptre_input
A=../../results/2026-04-14_refine_element_boundaries/elements_bed_final.bed

cut -f1-4 "${A}" > ${OUTDIR}/elements_bed_final_4col.bed

bedtools intersect -a "${OUTDIR}/elements_bed_final_4col.bed" -b ../../annotations/genes/gencode.v43.pc.genes.bed -wo -nonamecheck > ${OUTDIR}/el_vs_gene.wo.tsv
bedtools intersect -a "${OUTDIR}/elements_bed_final_4col.bed" -b ../../annotations/genes/gencode.v43.pc.exons.bed -wo -nonamecheck > ${OUTDIR}/el_vs_exon.wo.tsv
bedtools intersect -a "${OUTDIR}/elements_bed_final_4col.bed" -b ../../annotations/genes/gencode.v43.pc.cds.bed  -wo -nonamecheck > ${OUTDIR}/el_vs_cds.wo.tsv
bedtools intersect -a "${OUTDIR}/elements_bed_final_4col.bed" -b ../../annotations/genes/gencode.v43.pc.utr.bed  -wo -nonamecheck > ${OUTDIR}/el_vs_utr.wo.tsv
