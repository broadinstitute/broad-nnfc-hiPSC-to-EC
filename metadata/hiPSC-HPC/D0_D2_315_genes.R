###########################################################
#PRIMER DESIGN FOR D0-D2 SCREEN TARGET GENES WITH CONTROLS
###########################################################

#Assign the paths where primer3 (use ver 2.5 or higher) and blast tools (blastn, makeblastdb) are saved 
Sys.setenv(PATH = paste("/Users/jray/primer3-2.5.0/src", Sys.getenv("PATH"), 
                        sep=":"))

Sys.setenv(PATH = paste("/Users/jray/ncbi-blast-2.14.0+/bin", Sys.getenv("PATH"), 
                        sep = ":"))

#load the required libraries
library(TAPseq)
library(GenomicRanges)
library(BiocParallel)

#Make sure TAPseq is calling the primer3 and blastn from the correct paths

#Download gencode.v32.annotation.gtf from UCSC (hg38)
#convert it to a GRanges dataset
all_chr_genes <- rtracklayer::import('/Users/jray/TAP-seq/gencode.v32.annotation.gtf')

#head(all_chr_genes)

library(readxl)
#Remove any bad transcripts. Make a list of bad transcripts of the desired genes with isoforms not matching the refseq track
rm_txs <- read_xlsx("/Users/jray/TAP-seq/iPSC-HPC_D0_D2_screen/Bad_transcripts.xlsx")
# get the ENST numbers
rm_txs_genes <- rm_txs$transcript_id
rm_txs_genes
rm_txs_genes <- rm_txs_genes[!is.na(rm_txs_genes)]
rm_txs_genes
# remove those isoforms
all_chr_genes <- all_chr_genes[!all_chr_genes$transcript_id %in% rm_txs_genes]

# Check if any transcript_id from rm_txs_genes is still in all_chr_genes
#remaining_genes <- all_chr_genes$transcript_id[all_chr_genes$transcript_id %in% rm_txs_genes]

# If remaining_genes is empty, the removal was successful
#if (length(remaining_genes) == 0) {
#  print("All specified transcript_ids have been successfully removed.")
#} else {
#  print("Some transcript_ids were not removed:")
# print(remaining_genes)
#}




#Extract exons only
all_chr_genes_exons <- all_chr_genes[which(all_chr_genes$type == 'exon')]


# convert to GRangesList containing annotated exons per gene.
d0_target_genes <- split(all_chr_genes_exons, f = all_chr_genes_exons$gene_name)
d2_target_genes <- split(all_chr_genes_exons, f = all_chr_genes_exons$gene_name)

#extract the annotated target genes names for d0+d2
d0_genes <- read.csv("/Users/jray/TAP-seq/iPSC-HPC_D0_D2_screen/For_315_genes/d0_d2_with_old_ctrl_genes_for_design.csv")
d0_gene_list <- d0_genes$gene
d0_gene_list
d0_target_genes <- d0_target_genes[d0_gene_list]
d0_target_genes

#extract the annotated target genes names for d2 only
d2_genes <- read.csv("/Users/jray/TAP-seq/iPSC-HPC_D0_D2_screen/d2_only_genes_for_primer_design.csv")
d2_gene_list <- d2_genes$gene
d2_gene_list
d2_target_genes <- d2_target_genes[d2_gene_list]
d2_target_genes

#write.csv(target_genes, "/Users/jray/TAP-seq/hWAT-TF-TAP-primer-design/selec_target_genes.csv" )

#load the D0 WTC11 100TF perturb seq data (source: gs://landerlab-20210210-elisa-ipsc-villages/2022_06_25_WTC_100TFs_10x/results/D0_WTC_100TFs/possorted_genome_bam.bam)
d0_bam <- "/Users/jray/TAP-seq/iPSC-HPC_D0_D2_screen/bam_files/D0_100TF_WTC11_perturbseq/2022_06_25_WTC_100TFs_10x_results_D0_WTC_100TFs_possorted_genome_bam.bam"

#load the D2 WTC11 100TF perturb seq data (source: gs://landerlab-20210210-elisa-ipsc-villages/2022_06_25_WTC_100TFs_10x/results/D2_WTC_100TFs/possorted_genome_bam.bam)
d2_bam <- "/Users/jray/TAP-seq/iPSC-HPC_D0_D2_screen/bam_files/D2_100TF_WTC11_perturbseq/2022_06_25_WTC_100TFs_10x_results_D2_WTC_100TFs_possorted_genome_bam.bam"

#Check chr names in the bam file . This could be a. cause of error in running  infer_polyA_sites
#library(Rsamtools)
#test_bam <- scanBamHeader(d0_bam)
#bam_chr_names <- names(test_bam[[1]]$targets)
#print(bam_chr_names)

#If they do not have "chr" in their names, then instead of changing the bam, change the target_genes object
# Remove 'chr' from the chromosome names in target_genes
seqlevels(d0_target_genes) <- gsub("^chr", "", seqlevels(d0_target_genes))
d0_target_genes

seqlevels(d2_target_genes) <- gsub("^chr", "", seqlevels(d2_target_genes))
d2_target_genes


# register backend for parallelization
register(MulticoreParam(workers = 5))

# infer polyA sites for d0 genes from 10x data
d0_polyA_sites <- inferPolyASites(d0_target_genes, bam = d0_bam, polyA_downstream = 100,
                               wdsize = 200, min_cvrg = 1, parallel = TRUE)

#Manually inspect polyA site predictions and remove any obvious  false positives
#split polyA sites by gene (if you want to edit pA sites)
polyA_sites_split <- split(d0_polyA_sites, f = names(d0_polyA_sites))


#remove polyA sites for a gene. For the gene with removed pA sites, the 3' end of the last exon is used to truncate the transcript.
polyA_sites_split <- polyA_sites_split[! names(polyA_sites_split) %in% c("AP002495.1", "CADM2", "ESCO1", "HLA-E", "LRRC4C", "LSAMP", "RBMS3", "TTC25", "ZCWPW2", "ZNF608", "GRIA4", "PEX5L")]

#remove or keep specific polyA sites for a gene. Numbering starts from left to right
polyA_sites_split$ATN1 <- polyA_sites_split$ATN1[-c(2,3)]
polyA_sites_split$AUTS2 <- polyA_sites_split$AUTS2[15]
polyA_sites_split$B3GNTL1 <- polyA_sites_split$B3GNTL1[2]
polyA_sites_split$BECN1 <- polyA_sites_split$BECN1[1]
polyA_sites_split$C1orf50 <- polyA_sites_split$C1orf50[1]
polyA_sites_split$CDC20 <- polyA_sites_split$CDC20[2]
polyA_sites_split$CNP <- polyA_sites_split$CNP[2]
polyA_sites_split$CTNND1 <- polyA_sites_split$CTNND1[5]
polyA_sites_split$DCXR <- polyA_sites_split$DCXR[1]
polyA_sites_split$EMG1 <- polyA_sites_split$EMG1[2]
polyA_sites_split$FOXA2 <- polyA_sites_split$FOXA2[1]
polyA_sites_split$GNL1 <- polyA_sites_split$GNL1[2]
polyA_sites_split$KNOP1 <- polyA_sites_split$KNOP1[1]
polyA_sites_split$LRRC8B <- polyA_sites_split$LRRC8B[1]
polyA_sites_split$MLX <- polyA_sites_split$MLX[1]
polyA_sites_split$MRPL43 <- polyA_sites_split$MRPL43[2]
polyA_sites_split$NANOG <- polyA_sites_split$NANOG[1]
polyA_sites_split$NARF <- polyA_sites_split$NARF[-c(7)]
polyA_sites_split$OGFOD3<- polyA_sites_split$OGFOD3[5]
polyA_sites_split$PCYT2 <- polyA_sites_split$PCYT2[-c(1,2)]
polyA_sites_split$RIMKLA <- polyA_sites_split$RIMKLA[2]
polyA_sites_split$RIMKLB <- polyA_sites_split$RIMKLB[3]
polyA_sites_split$RNF121 <- polyA_sites_split$RNF121[1]
polyA_sites_split$SAR1A <- polyA_sites_split$SAR1A[3]
polyA_sites_split$SGPL1 <- polyA_sites_split$SGPL1[2]
polyA_sites_split$SHANK2 <- polyA_sites_split$SHANK2[2]
polyA_sites_split$SLC16A3 <- polyA_sites_split$SLC16A3[1]
polyA_sites_split$SZT2 <- polyA_sites_split$SZT2[11]
polyA_sites_split$TBCD <- polyA_sites_split$TBCD[2]
polyA_sites_split$TYSND1 <- polyA_sites_split$TYSND1[3]
polyA_sites_split$UNC5B <- polyA_sites_split$UNC5B[3]
polyA_sites_split$VARS <- polyA_sites_split$VARS[2]
polyA_sites_split$TSHZ3 <- polyA_sites_split$TSHZ3[1]

#check the final edited list to make sure the desired edits are there
#View(as.data.frame(polyA_sites_split))

#revert the pruned pA sites to GRanges object
d0_polyA_sites <- unlist(polyA_sites_split)

# infer polyA sites for d2 genes from 10x data
d2_polyA_sites <- inferPolyASites(d2_target_genes, bam = d2_bam, polyA_downstream = 100,
                                  wdsize = 200, min_cvrg = 1, parallel = TRUE)


#seqlevels(d0_polyA_sites)
#seqlevels(d2_polyA_sites)

#split d2 polyA sites by gene (if you want to edit pA sites)
d2_polyA_sites_split <- split(d2_polyA_sites, f = names(d2_polyA_sites))

#remove polyA sites for a gene. For the gene with removed pA sites, the 3' end of the last exon is used to truncate the transcript.
d2_polyA_sites_split <- d2_polyA_sites_split[! names(d2_polyA_sites_split) %in% c("HAS2")]

#remove or keep specific polyA sites for a gene. Numbering starts from left to right
d2_polyA_sites_split$EGFLAM <- d2_polyA_sites_split$EGFLAM[1]
d2_polyA_sites_split$PLXNA2 <- d2_polyA_sites_split$PLXNA2[4]

#check the final edited list to make sure the desired edits are there
#View(as.data.frame(d2_polyA_sites_split))

#revert the pruned pA sites to GRanges object
d2_polyA_sites <- unlist(d2_polyA_sites_split)

#Combine the polyA sites and target genes GRanges objects of d0 and d2
polyA_sites <- c(d0_polyA_sites, d2_polyA_sites)
target_genes <- c(d0_target_genes, d2_target_genes)


#seqlevels(polyA_sites)
#seqlevels(target_genes)


library(rtracklayer)
export(polyA_sites, con = "/Users/jray/TAP-seq/iPSC-HPC_D0_D2_screen/outputs/run_6/315_final_genes/final_315_edit_polyA_sites.bed", format = "bed")
export(unlist(target_genes), con = "/Users/jray/TAP-seq/iPSC-HPC_D0_D2_screen/outputs/run_6/315_final_genes/final_315_edit_target_genes.gtf", format = "gtf")

# truncate transcripts at inferred polyA sites
truncated_txs <- truncateTxsPolyA(target_genes, polyA_sites = polyA_sites, parallel = TRUE)
truncated_txs
export(unlist(truncated_txs), con = "/Users/jray/TAP-seq/iPSC-HPC_D0_D2_screen/outputs/run_6/315_final_genes/final_315_edit_truncated_txs.gtf", format = "gtf")


library(BSgenome)

# human genome BSgenome object (needs to be installed from Bioconductor)
#library(BiocManager)
#install("BSgenome.Hsapiens.UCSC.hg19")
hg38 <- getBSgenome("BSgenome.Hsapiens.UCSC.hg38")

truncated_txs

# Add "chr" prefix to the chromosome names in truncated_txs
seqlevels(truncated_txs) <- paste0("chr", seqlevels(truncated_txs))
truncated_txs

# get sequence for all truncated transcripts
txs_seqs <- getTxsSeq(truncated_txs, genome = hg38)


#DESIGN OUTER PRIMERS
# create TAPseq IO for outer forward primers from truncated transcript sequences
outer_primers <- TAPseqInput(txs_seqs, target_annot = truncated_txs,
                             product_size_range = c(350, 500), primer_num_return = 5)
# design 5 outer primers for each target gene
outer_primers <- designPrimers(outer_primers)

#outer_primers
# check outer_primers
#list_op <- tapseq_primers(outer_primers)
#list_op
#View(as.data.frame(list_op))

# get pcr products for specific genes
#pcr_products(outer_primers$GATA1)

#Blast primers against whole genome and transcriptome. Build a blast database the first time. 
#The database can be saved in a location and used for multiple primer designs whithout having to rebuild it everytime.
library(BSgenome)

# human genome BSgenome object (needs to be installed from Bioconductor)
hg38 <- getBSgenome("BSgenome.Hsapiens.UCSC.hg38")


# download and import gencode hg38 annotations
url <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.annotation.gtf.gz"
annot <- import(url, format = "gtf")
#If there is a download error then use the following
#annot <- import('/Users/jray/TAP-seq/gencode.v32.annotation.gtf', format = "gtf")

# extract exon annotations for protein-coding genes to build transcripts
tx_exons <- annot[annot$type == "exon" & annot$gene_type == "protein_coding"]

#Make sure the correct path to makeblastdb is assigned
#options(TAPseq.makeblastdb = "/Users/jray/ncbi-blast-2.14.0+/bin/makeblastdb")

# create blast database
blastdb <- file.path(tempdir(), "blastdb") 
#Can run print(blastdb) to check the path and access it by running "open /path/to/blastdb" in terminal
createBLASTDb(genome = hg38, annot = tx_exons, blastdb = blastdb)

# now we can blast our outer primers against the created database (run this for small gene panels)
#outer_primers <- blastPrimers(outer_primers, blastdb = blastdb, max_mismatch = 0, min_aligned = 0.75)

#Run the blast step in sequential batches (for large gene panels)
# 1. Split outer_primers into batches of a maximum of 10 genes per batch
batch_size <- 10  # Adjust batch size as needed
primer_batches <- split(outer_primers, ceiling(seq_along(outer_primers) / batch_size))

# Verify the batches
lapply(primer_batches, function(batch) class(batch))

#Verfy the gene list in each batch
# Iterate over each batch in primer_batches
#for (i in seq_along(primer_batches)) {
  # Print batch number
#  cat("Batch", i, "\n")
  
  # Extract the current batch
#  current_batch <- primer_batches[[i]]
  
  # Iterate over each TsIO object in the batch
#  for (j in seq_along(current_batch)) {
    # Extract the TsIO object (which contains the gene name)
#    tsio_object <- current_batch[[j]]
    
    # Extract and print the gene name (replace with the correct accessor if needed)
#    gene_name <- tsio_object@sequence_id  # Assuming sequence_id contains gene names
#    cat("Gene:", gene_name, "\n")
#  }
  
  # Add a line break between batches for readability
#  cat("\n")
#}


# 2. Function to run blastPrimers on each batch and time it
blast_results_batches <- lapply(seq_along(primer_batches), function(batch_num) {
  cat("Running BLAST for batch", batch_num, "of", length(primer_batches), "...\n")
  
  # Run blastPrimers on this batch
  result <- tryCatch({
    blastPrimers(primer_batches[[batch_num]], blastdb = blastdb, max_mismatch = 0, min_aligned = 0.75)
  }, error = function(e) {
    cat("Error in batch", batch_num, ":", e$message, "\n")
    NULL  # Return NULL if there's an error
  })
  
  cat("Batch", batch_num, "completed.\n")
  
  return(result)
})

# 3. Combine all valid TsIOList objects from the batch results
valid_batches <- Filter(function(batch) !is.null(batch) && inherits(batch, "TsIOList"), blast_results_batches)

if (length(valid_batches) > 0) {
  combined_results <- do.call(c, valid_batches)
  cat("Successfully combined all valid batches.\n")
  
  # 4. Save the combined results as outer_primers
  outer_primers <- combined_results
  cat("Saved combined results as outer_primers.\n")
  
} else {
  cat("No valid batches to combine.\n")
}

outer_primers
  
# the primers now contain the number of estimated off-targets
#list_op_off <- tapseq_primers(outer_primers)


#View(as.data.frame(list_op_off))

# select best primer per target gene based on the number of potential off-targets
best_outer_primers <- pickPrimers(outer_primers, n = 1, by = "off_targets")

# each object now only contains the best primer
#list_op_best <- tapseq_primers(best_outer_primers)
################################################################################
#DESIGN INNER PRIMERS
# create new TsIO objects for inner primers, note the different product size
inner_primers <- TAPseqInput(txs_seqs, target_annot = truncated_txs,
                             product_size_range = c(150, 300), primer_num_return = 5)


# design inner primers
inner_primers <- designPrimers(inner_primers)
# blast inner primers (run this for small gene panels)
#inner_primers <- blastPrimers(inner_primers, blastdb = blastdb, max_mismatch = 0, min_aligned = 0.75)

#Run the blast step in sequential batches (for large gene panels) 
# 1. Split inner_primers into batches of a maximum of 10 genes per batch
batch_size <- 10  # Adjust the batch size if needed
ip_primer_batches <- split(inner_primers, ceiling(seq_along(inner_primers) / batch_size))

# 2. Function to run blastPrimers on each batch and time it
ip_blast_results_batches <- lapply(seq_along(ip_primer_batches), function(batch_num) {
  cat("Running BLAST for batch", batch_num, "of", length(ip_primer_batches), "...\n")
  
  # Run blastPrimers on this batch
  result <- tryCatch({
    blastPrimers(ip_primer_batches[[batch_num]], blastdb = blastdb, max_mismatch = 0, min_aligned = 0.75)
  }, error = function(e) {
    cat("Error in batch", batch_num, ":", e$message, "\n")
    NULL  # Return NULL if there's an error
  })
  
  cat("Batch", batch_num, "completed.\n")
  
  return(result)
})

# 3. Combine all valid TsIOList objects from the batch results
ip_valid_batches <- Filter(function(batch) !is.null(batch) && inherits(batch, "TsIOList"), ip_blast_results_batches)

if (length(ip_valid_batches) > 0) {
  ip_combined_results <- do.call(c, ip_valid_batches)
  cat("Successfully combined all ip valid batches.\n")
  
  # 4. Save the combined results as inner_primers
  inner_primers <- ip_combined_results
  cat("Saved combined results as inner_primers.\n")
  
} else {
  cat("No valid batches to combine.\n")
}
inner_primers
# the primers now contain the number of estimated off-targets
#list_ip_off <- tapseq_primers(inner_primers)

# pick best primer per target gene
best_inner_primers <- pickPrimers(inner_primers, n = 1, by = "off_targets")
# each object now only contains the best primer
#list_ip_best <- tapseq_primers(best_inner_primers)
#save in a csv file
#write.csv(list_ip_best,"/Users/jray/TAP-seq/Primer_list/Inner_all_best.csv")
#################################################################################
#CHECK MULTIPLEX COMPATIBILTY
# use checkPrimers to run Primer3's "check_primers" functionality for every possible primer pair
outer_comp <- checkPrimers(best_outer_primers)
inner_comp <- checkPrimers(best_inner_primers)

library(dplyr)
#install.packages('ggplot2')
library(ggplot2)

# merge outer and inner complementarity scores into one data.frame
comp <- bind_rows(outer = outer_comp, inner = inner_comp, .id = "set")

# add variable for pairs with any complemetarity score higher than 47
comp <- comp %>%
  mutate(high_compl = if_else(primer_pair_compl_any_th > 47 | primer_pair_compl_end_th > 47,
                              true = "high", false = "ok")) %>% 
  mutate(high_compl = factor(high_compl, levels = c("ok", "high")))

# plot complementarity scores (make an interactive plot)
p <- ggplot(comp, aes(x = primer_pair_compl_any_th, y = primer_pair_compl_end_th)) +
  facet_wrap(~set, ncol = 2) +
  geom_point(aes(color = high_compl), alpha = 0.25) +
  scale_color_manual(values = c("black", "red"), drop = FALSE) +
  geom_hline(aes(yintercept = 47), colour = "darkgray", linetype = "dashed") +
  geom_vline(aes(xintercept = 47), colour = "darkgray", linetype = "dashed") +
  labs(title = "Complementarity scores TAP-seq primer combinations",
       color = "Complementarity") +
  theme_bw()

library(plotly)
ggplotly(p)

#View the compatibility table to identify problematic pairs
#View(comp)
#write.csv(comp,"/Users/jray/TAP-seq/iPSC-HPC_D0_D2_screen/outputs/run_4/all_edit_genes_outs/all_edit_compatibility.csv")

#EXPORT PRIMERS
# create data.frames for outer and inner primers
outer_primers_df <- primerDataFrame(best_outer_primers)
inner_primers_df <- primerDataFrame(best_inner_primers)

# the resulting data.frames contain all relevant primer data
#View(outer_primers_df)
#View(inner_primers_df)
#save in a csv file
write.csv(outer_primers_df,"/Users/jray/TAP-seq/iPSC-HPC_D0_D2_screen/outputs/run_6/315_final_genes/final_315_edit_Outer.csv")
write.csv(inner_primers_df,"/Users/jray/TAP-seq/iPSC-HPC_D0_D2_screen/outputs/run_6/315_final_genes/final_315_edit_Inner.csv")

# create BED tracks for outer and inner primers with custom colors
outer_primers_track <- createPrimerTrack(best_outer_primers, color = "steelblue3")
inner_primers_track <- createPrimerTrack(best_inner_primers, color = "goldenrod1")

# the output data.frames contain lines in BED format for every primer
#head(outer_primers_track)
#head(inner_primers_track)
# export tracks to .bed files ("" writes to the standard output, replace by a file name)
exportPrimerTrack(outer_primers_track, con = "/Users/jray/TAP-seq/iPSC-HPC_D0_D2_screen/outputs/run_6/315_final_genes/final_315_edit_Outer.bed")
exportPrimerTrack(inner_primers_track, con = "/Users/jray/TAP-seq/iPSC-HPC_D0_D2_screen/outputs/run_6/315_final_genes/final_315_edit_Inner.bed")
#Workflow complete#
