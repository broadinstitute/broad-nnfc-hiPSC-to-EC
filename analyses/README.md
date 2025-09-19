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

- **[positive_pairs_hg38_readout_genes.tsv](../../results/2025-09-15/positive_pairs_hg38_readout_genes.tsv)**  
  TSV file containing positive element-gene pairs.  
  **Columns:**  
  - `element_name`: Name of the element  
  - `ensembl_id`: Ensembl ID of the TSS being targeted

- **[discovery_pairs_no_pos_hg38_readout_genes.tsv](../../results/2025-09-15/discovery_pairs_no_pos_hg38_readout_genes.tsv)**  
  TSV file with discovery pairs (excluding positive pairs and TSS overlaps).  
  **Columns:**  
  - `element_name`: Name of the element  
  - `ensembl_id`: Ensembl ID to test for interaction

- **[grna_non_targeting_input.tsv](../../results/2025-09-15/grna_non_targeting_input.tsv)**  
  TSV file listing non-targeting gRNAs used in the experiment.  
  **Columns:**  
  - `grna_id`: gRNA identifier  
  - `label`: "non-targeting" (fixed value)

- **[d0_seurat_object_filtered_norm.rds](../../results/2025-09-18/d0_seurat_object_filtered_norm.rds)**  
  RDS file containing a Seurat object with filtered and normalized data for D0 hiPSC.

- **[mean_cpm_per_gene_filtered_cells.tsv](../../results/2025-09-18/mean_cpm_per_gene_filtered_cells.tsv)**  
  TSV file with mean CPM per gene across all filtered cells (lenient filtering).  
  **Columns:**  
  - `gene_symbol`: Gene symbol  
  - `ensembl_id`: Ensembl gene ID  
  - `mean_cpm`: Mean CPM value across all filtered cells  
  - `is_readout_gene`: Indicates if the gene is a readout gene (TRUE/FALSE)

- **[high_cpm_non_readout_genes.tsv](../../results/2025-09-18/high_cpm_non_readout_genes.tsv)**  
  TSV file listing non-readout genes with high expression (mean CPM > 100).  
  **Columns:**  
  - `gene_symbol`: Gene symbol  
  - `ensembl_id`: Ensembl gene ID  
  - `mean_cpm`: Mean CPM value

- **[gene_annotations.csv](../../results/2025-09-18_gene_annotations/gene_annotations.csv)**  
  CSV file with annotations for the high CPM non-readout genes. 
  **Columns:**  
  - `gene_symbol`: Gene symbol  
  - `name`: This is the full gene name
  - `summary`: Description of the function
  
## hiPSC → EC differentiation large screen

### Code
- [2025-09-10](./2025-09-10/): First look at the hiPSC-EC large screen

### Results