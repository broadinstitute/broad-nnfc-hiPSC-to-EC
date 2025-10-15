# Checking inner primer usage

In the DC-TAP readout we see a lot of UMIs coming from genes that are not the intended target of the primers. This analysis checks how many UMIs come from transcripts that are actually targeted by the inner primers.

The script `check_inner_primer_usage.py` takes as input the BAM file with the aligned reads and a CSV file with the intended target genes and their inner primers. It checks for each read if it is assigned to an intended target gene and what inner primer is found in the read.

## Algorithm Efficiency

The script uses an optimized dictionary lookup approach:

- For each read, it makes at most `(max_primer_length - min_primer_length + 1)` O(1) dictionary queries
- This is much more efficient than iterating through all primers for each read
- Handles variable primer lengths automatically by testing sequence slices of different lengths

## Input files

- `possorted_genome_bam.bam`: BAM file with the aligned reads (output of Cell Ranger)
- `[metadata_genes_final.csv](../../metadata/metadata_genes_final.csv)`: CSV file with the intended target genes (column id: `intended_target_gene_id`) and their inner primers (column id: `inner_primer`) (output of `analyses/00-2025-10-08_create_gene_metadata/README.md`)

## Usage

```bash
# Basic usage - search for primers at 5' end
python check_inner_primer_usage.py possorted_genome_bam.bam metadata_genes_final.csv

# Search for primers at 3' end
python check_inner_primer_usage.py possorted_genome_bam.bam metadata_genes_final.csv --3-prime

# Reverse complement primers before searching
python check_inner_primer_usage.py possorted_genome_bam.bam metadata_genes_final.csv --reverse-complement

# Analyze primer positions and orientations
python check_inner_primer_usage.py possorted_genome_bam.bam metadata_genes_final.csv --analyze-positions
```

## Output files

- `primer_usage_report.txt`: Summary statistics and top gene assignments
- `primer_usage_report_detailed.tsv`: Detailed results for each read with primer found
- `primer_position_analysis.txt`: Analysis of primer positions and orientations (when using `--analyze-positions`)

