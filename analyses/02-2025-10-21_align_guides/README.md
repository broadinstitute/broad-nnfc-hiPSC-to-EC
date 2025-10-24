# Aligning CRISPR guides to the genome

## Generating FASTQ files from guide sequences
The initial exploration and the creation of the FASTQs files  was done using the marimo notebook [./create_guide_fastq.py](./create_guide_fastq.py).

The notebook generates three FASTQ files:
 - `ec_differentiation.fastq`: contains the guide sequences without PAM or leading G
 - `ec_differentiation_with_pam.fastq`: contains the guide sequences with PAM but without leading G
 - `ec_differentiation_with_pam_and_leading_g.fastq`: contains the guide sequences with PAM and leading G

## Aligning FASTQ files to the genome
The genome FASTA files were downloaded from ENCODE:
 - [ENCODE hg38](https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz)

Reads are aligned to the genome using the GEM3 mapper.

```bash
  ../../gem3-mapper/bin/gem-mapper \
  -I ../../annotations/ENCODE/hg38/gem_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.gem \
  -i ec_differentiation_with_pam_and_leading_g.fastq \
  -o ec_differentiation_with_pam_and_leading_g_hg38.sam \
  --sam-compact=false \
  --mapping-mode sensitive \
  --threads 8
```

## Filtering alignments and annotating guides

The script [filter_alignments_with_logs.py](filter_alignments_with_logs.py) is used to select valid alignments.
An alignment is considered valid if it meets the following criteria:
 - Maps to the primary assembly
 - NGG PAM must be aligned. Mismatches, insertions, deletions, and soft-clipping are not allowed. The mismatch for the 'N' nucleotide of the 'NGG' PAM is accepted.
 - The spacer must map without insertions or deletions
 - The leading 'G' can be softclipped

If multiple valid alignments for the same guide exist, all are kept.


```bash
python3 ../../analyses/02-2025-10-21_align_guides/filter_alignments_with_logs.py ec_differentiation_with_pam_and_leading_g_hg38_namesort.sam \
--require-pam-clean \
--emit-protospacer --drop-leading-g \
--out-bed ec_differentiation_with_pam_and_leading_g_hg38_filtered.bed \
--out-discarded discarded.tsv \
--out-unmapped unmapped.tsv \
--out-mismatch-csv mismatches.csv \
--plot-html mismatch_histogram.html
```

## Output Files

### BED file (`ec_differentiation_with_pam_and_leading_g_hg38_filtered.bed`)
Contains valid alignments in BED format with additional columns:

| Column | Description |
|--------|-------------|
| 1. chrom | Chromosome name |
| 2. chromStart | Start position (0-based) |
| 3. chromEnd | End position (0-based, exclusive) |
| 4. name | Guide read name |
| 5. score | Mapping quality (0-255) |
| 6. strand | Strand (+ or -) |
| 7. NM | Number of mismatches |
| 8. AS | Alignment score |

### Discarded alignments (`discarded.tsv`)
Contains mapped alignments that failed filtering criteria:

| Column | Description |
|--------|-------------|
| read_name | Guide read name |
| chromosome | Chromosome name |
| pos1 | Start position (1-based) |
| flag | SAM flag |
| mapq | Mapping quality |
| strand | Strand (+ or -) |
| cigar | CIGAR string |
| NM | Number of mismatches |
| AS | Alignment score |
| MD | MD tag (mismatch/deletion string) |
| reason | Reason for discarding (see below) |

**Discard reasons:**
- `discarded_tail_unaligned`: PAM region not properly aligned
- `discarded_tail_mismatch`: PAM region contains mismatches
- `discarded_protospacer_indel`: Protospacer region contains insertions or deletions

### Unmapped reads (`unmapped.tsv`)
Contains reads that did not map to the genome:

| Column | Description |
|--------|-------------|
| read_name | Guide read name |
| flag | SAM flag |

### Mismatch details (`mismatches.csv`)
Contains detailed mismatch information for analysis:

| Column | Description |
|--------|-------------|
| read_name | Guide read name |
| chromosome | Chromosome name |
| strand | Strand (+ or -) |
| read_length | Total query sequence length |
| query_index | Position of mismatch in query (0-based) |
| distance_from_pam | Distance from PAM end (0 = PAM adjacent) |
| nm | Total number of mismatches in alignment |
| as | Alignment score |
| reason | Whether alignment was kept or discarded |

### Mismatch histogram (`mismatch_histogram.html`)
Interactive HTML visualization showing mismatch frequency by position relative to PAM, faceted by alignment fate (kept vs discarded reasons).

