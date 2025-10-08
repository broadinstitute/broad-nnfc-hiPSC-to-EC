# Align the guides with PAM to the genome

To check wheter the guides map to expected locations they were mapped to the genome using the [gem aligner](https://github.com/smarco/gem3-mapper).

## Creating input FASTA and FASTQ files

First the FASTA and FASTQ files were created from the [220308_K562_Random_Screen_Crop.design](../../data/Design_documents/220308_K562_Random_Screen_Crop.design.txt) and [220308_WTC11_Random_Screen_Crop.design](../../data/Design_documents/220308_WTC11_Random_Screen_Crop.design.txt) files.

The script [fasta2fastq.awk](../../src/awk/fasta2fastq.awk) converts the format from FASTA to FASTQ. Max Phred33 score( `I` ) is used except for `N` nucleotides that get the lowest Phred33 score ( `!` ).

### Convert K562

```
grep K562 ../../legacy/Supplementary_tables/Table_S5_Guide_list.tsv | grep -v negative_control | grep -v safe_targeting | grep -wv DE | grep -wv tss_pos | cut -f4 > K562_names_guides_targeting_loci.txt

cut -f4,13 ../../data/220308_guide_design/220308_K562_Random_Screen_Crop.design.txt | awk 'NR>1{print ">"$1"\n"$2}' > K562_guides_with_PAM.fasta


awk -f ../../src/awk/fasta2fastq.awk K562_guides_with_PAM.fasta > K562_guides_with_PAM.fastq

bgzip K562_guides_with_PAM.fast*
```

### Convert WTC11

```
grep WTC11 ../../legacy/Supplementary_tables/Table_S5_Guide_list.tsv | grep -v negative_control | grep -v safe_targeting | grep -wv DE | grep -wv tss_pos | cut -f4 > WTC11_names_guides_targeting_loci.txt

cut -f4,13 ../../data/220308_guide_design/220308_WTC11_Random_Screen_Crop.design.txt | awk 'NR>1{print ">"$1"\n"$2}' > WTC11_guides_with_PAM.fasta

awk -f ../../src/awk/fasta2fastq.awk WTC11_guides_with_PAM.fasta > WTC11_guides_with_PAM.fastq

bgzip WTC11_guides_with_PAM.fast*
```

## Aligning FASTQ files to the genome

The genome FASTA files were downloaded from GENCODE:
 - [ENCODE hg38](https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz)
 - [ENCODE hg19 female](https://www.encodeproject.org/files/female.hg19/@@download/female.hg19.fasta.gz)

**Creating indexes for the genomes:**

*ENCODE hg38*

`../../../gem3-mapper/bin/gem-indexer -i genome/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz -o gem_index/GRCh38_no_alt_analysis_set_GCA_000001405.15 -t 8`

*ENCODE hg19 female*

`../../../gem3-mapper/bin/gem-indexer -i genome/female.hg19.fasta.gz -o gem_index/ENCODE_female.hg19 -t 8`

**Aligning guides to the genomes**

*GRCh38.p13*

```
../gem3-mapper/bin/gem-mapper \
  -I ../annotations/ENCODE/hg38/gem_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.gem \
  -i K562_guides_with_PAM.fastq.gz \
  -o K562_guides_with_PAM_hg38.sam \
  --mapping-mode sensitive \
  --threads 8

../gem3-mapper/bin/gem-mapper \
  -I ../annotations/ENCODE/hg38/gem_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.gem \
  -i WTC11_guides_with_PAM.fastq.gz \
  -o WTC11_guides_with_PAM_hg38.sam \
  --mapping-mode sensitive \
  --threads 8
```


*GRCh37.p13*

```
../gem3-mapper/bin/gem-mapper \
  -I ../annotations/ENCODE/hg19/gem_index/ENCODE_female.hg19.gem \
  -i K562_guides_with_PAM.fastq.gz \
  -o K562_guides_with_PAM_hg19.sam \
  --mapping-mode sensitive \
  --threads 8

../gem3-mapper/bin/gem-mapper \
  -I ../annotations/ENCODE/hg19/gem_index/ENCODE_female.hg19.gem \
  -i WTC11_guides_with_PAM.fastq.gz \
  -o WTC11_guides_with_PAM_hg19.sam \
  --mapping-mode sensitive \
  --threads 8
```



