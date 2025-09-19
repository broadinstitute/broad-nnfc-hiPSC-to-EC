# Index
- [hiPSC → HPC differentiation large screen](#hipsc--hpc-differentiation-large-screen)
- [hiPSC → EC differentiation large screen](#hipsc--ec-differentiation-large-screen)


## hiPSC → HPC differentiation large screen

### Code

- [2025-09-02](./2025-09-02/): DEPRECATED: Preliminary analysis with inaccurate pairs. ~~First look at the D0 iPSC -> HPC large screen.~~
- [2025-09-15](./2025-09-15/): Create discovery and positive pairs
- [2025-09-18](./2025-09-18/): Load data for D0 hiPSC -> HPC and preprocess
- [2025-09-18_gene_annotations](./2025-09-18_gene_annotations/): WARNING: This is an alpha version. Given a list of genes, get their annotations (biotype, entrez id, etc) using my_gene python library.
- [2025-09-19](./2025-09-19/): SCEPTRE analysis for D0 hiPSC -> HPC

### Results
- [positive_pairs_hg38_readout_genes.tsv](../../results/2025-09-15/positive_pairs_hg38_readout_genes.tsv): TSV file, first column is element name, ensembl id of the tss being targeted.
- [discovery_pairs_no_pos_hg38_readout_genes.tsv](../../results/2025-09-15/discovery_pairs_no_pos_hg38_readout_genes.tsv): TSV file, first column is element name, ensembl id we want to test for interaction. All the elemnts overlapping a TSS have been removed from this file.
- [grna_non_targeting_input.tsv](../../results/2025-09-15/grna_non_targeting_input.tsv): TSV file, list of non-targeting gRNAs used in the experiment.
- [d0_seurat_object_filtered_norm.rds](../../results/2025-09-18/d0_seurat_object_filtered_norm.rds): RDS file, Seurat object with filtered and normalized data for D0 hiPSC
- [mean_cpm_per_gene_filtered_cells.tsv](../../results/2025-09-18/mean_cpm_per_gene_filtered_cells.tsv): TSV file, mean CPM per gene across all filtered cells. Columns are gene_symbol, ensembl_id, mean_cpm, is_readout_gene. This is a targeted transcriptome experiment, the last column indicates if the gene is one of the readout genes.
- [high_cpm_non_readout_genes.tsv](../../results/2025-09-18/high_cpm_non_readout_genes.tsv): TSV file, genes that are not readout genes but have high expression (mean CPM > 100) in the filtered cells. Columns are gene_symbol, ensembl_id, mean_cpm.
- [gene_annotations.csv](../../results/2025-09-18_gene_annotations/gene_annotations.csv): CSV file, annotations for all the genes in the dataset (readout and non-readout). Columns are gene_symbol, ensembl_id, entrez_id, biotype, name.

## hiPSC → EC differentiation large screen

### Code
- [2025-09-10](./2025-09-10/): First look at the hiPSC-EC large screen

### Results