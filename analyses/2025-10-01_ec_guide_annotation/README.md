# Generating FASTQs
The initial exploration and the creation of the FASTQs files  was done using the marimo notebook [./create_guide_fastq.pynb](./create_guide_fastq.pynb).

The notebook generates three FASTQ files:
 - `ec_differentiation.fastq`: contains the guide sequences without PAM or leading G
 - `ec_differentiation_with_pam.fastq`: contains the guide sequences with PAM but without leading G
 - `ec_differentiation_with_pam_and_leading_g.fastq`: contains the guide sequences with PAM and leading G

# Aligning FASTQ files to the genome
The genome FASTA files were downloaded from ENCODE:
 - [ENCODE hg38](https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz)

and aligned using GEM3 mapper.

```bash
../../gem3-mapper/bin/gem-mapper \
  -I ../../annotations/ENCODE/hg38/gem_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.gem \
  -i ec_differentiation.fastq \
  -o ec_differentiation_hg38.sam \
  --mapping-mode sensitive \
  --threads 8

  ../../gem3-mapper/bin/gem-mapper \
  -I ../../annotations/ENCODE/hg38/gem_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.gem \
  -i ec_differentiation_with_pam.fastq \
  -o ec_differentiation_with_pam_hg38.sam \
  --mapping-mode sensitive \
  --threads 8

  ../../gem3-mapper/bin/gem-mapper \
  -I ../../annotations/ENCODE/hg38/gem_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.gem \
  -i ec_differentiation_with_pam_and_leading_g.fastq \
  -o ec_differentiation_with_pam_and_leading_g_hg38.sam \
  --mapping-mode sensitive \
  --threads 8
  ```

  # Converting SAM to bed
  [convert_sam_with_xa_to_bed.py](/Users/emattei/GitHub/broad-nnfc-dc-tap/src/python/convert_sam_with_xa_to_bed.py) was used to convert the SAM files to BED format, taking into account the XA tag for multiple alignments.

  This BED file contains the following columns:
  - chrom: chromosome
  - start: start position
  - end: end position
  - guide_id: guide identifier
  - mapq: mapping quality
  - strand: strand
  - NM: number of mismatches
  - AS: alignment score
  - MD: mismatch description
  - type: type of alignment (primary or alternative)