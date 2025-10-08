# Annotate gene primers

Check if primers match the intended target genes using pysam to read a FASTA file of transcripts.

### Step 1

Install the required packages:
```
pip install csv, pysam, tqdm
```

Run the script:
```
python3 annotate_gene_primers.py
```

The script will take a csv with the inner primers used in the experiment and find all the transcripts that match the primers. It will then output a csv with the matches found.

Output file:

[primer_matches_updated_names.csv](../../results/2025-10-08_create_gene_metadata/intermediates/primer_matches_updated_names.csv)

Output columns:

`outer_primer`: The outer primer sequence used in the experiment.
`inner_primer`: The inner primer sequence used in the experiment.
`outer_distance_from_transcript_end`: Distance of the outer primer from the end of the transcript.
`inner_distance_from_transcript_end`: Distance of the inner primer from the end of the transcript.
`distance_between_primers`: Distance between the outer and inner primers if both are found.
`primers_found`: One of "both", "outer", "inner" indicating which primers were found.
`is_intended_target`: Boolean indicating if the primers match the intended target gene.
`intended_target_gene_symbol`: The gene symbol of the intended target.
`transcript_id`: The transcript ID where the primers were found.
`gene_id`: The gene ID where the primers were found.
`gene_symbol`: The gene symbol where the primers were found.
`gene_type`: The gene type where the primers were found.

### Step 2

Use the create gene metadata script to create a metadata file for all genes in the transcriptome.

Output files:

#### Metadata gene file
[metadata_genes_final.csv](../../metadata/metadata_genes_final.csv)
Columns:
 - intended_target_gene_id: Ensembl gene ID
 - intended_target_gene_symbol: Gene symbol
 - inner_primer: Inner primer sequence
 - outer_primer: Outer primer sequence
 - inner_primer_adapter: Inner primer sequence with adapters
 - is_central_gene: Is central gene (True/False)
 - locus_central_gene: Locus central gene the gene belongs to (if applicable)
 - well: Well position they came from the commercial provider
 - batch: Batch information


[ec_screen_16_loci_2mb_info_final.tsv](../../metadata/ec_screen_16_loci_2mb_info_final.tsv)
[ec_screen_16_loci_2mb_locus_final.bed](../../metadata/ec_screen_16_loci_2mb_locus_final.bed)
