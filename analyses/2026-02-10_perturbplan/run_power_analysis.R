#!/usr/bin/env Rscript

## Compute statistical power to detect a specific effect size on expression for all
## perturbation-gene pairs to test using perturbplan

## Process input arguments
library(optparse)

# define command line arguments
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL, 
              help = "Path to input MuData file (mudata.h5mu)", metavar = "file"),
  make_option(c("-o", "--output"), type = "character", default = NULL, 
              help = "Path to output MuData file (mudata.h5mu)", metavar = "file"),
  make_option(c("-e", "--effect_size"), type = "double", default = NULL,
              help = "Fold change effect size on expression, for which power will be computed",
              metavar = "real"),
  make_option(c("-s", "--effect_size_sd"), type = "double", default = 0.13,
              help = "Effect size standard deviation to model guide-guide variability",
              metavar = "real"),
  make_option(c("-p", "--pvalue_cutoff"), type = "double", default = 0.001,
              help = "Nominal p-value cutoff used in power analysis", metavar = "real"),
  make_option(c("-t", "--threads"), type = "integer", default = 1, 
              help = "Number of CPU threads to use for parallelization", metavar = "integer"),
  make_option("--control_group", type = "character", default = "complement",
              help = paste("String specifying the control group to use in the differential",
                "expression analysis, either 'complement' or 'nt_cells'"), metavar = "character"),
  make_option("--side", type = "character", default = "both",
              help = "The sidedness of the test, one of 'left', 'right', or 'both'",
              metavar = "character"),
  make_option("--n_nonzero_trt_thresh", type = "integer", default = 7,
              help = "Number of nonzero treatment cells a pair must contain for it to be retained",
              metavar = "character"),
  make_option("--n_nonzero_cntrl_thresh", type = "integer", default = 7,
              help = "Number of nonzero control cells a pair must contain for it to be retained",
              metavar = "character")      

)

# parse arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# function to check for required arguments
check_required_args <- function(arg, opt, opt_parser) {
  if (is.null(opt[[arg]])) {
    print_help(opt_parser)
    stop(arg, " argument is required!", call. = FALSE)
  }
}

# check that all required parameters are provided
required_args <- c("input", "output", "effect_size")
for (arg in required_args) {
  check_required_args(arg, opt = opt, opt_parser = opt_parser)
}

# save.image("power_analysis.rda")
# stop()

## Define functions --------------------------------------------------------------------------------

# required packages
suppressPackageStartupMessages({
  library(MuData)
  library(MultiAssayExperiment)
  library(BiocParallel)
  library(perturbplan)
  library(sceptre)
  library(sceptreIGVF)
  library(data.table)
})

# run sceptre precomputation an all genes to get dispersion values
calculate_sceptre_dispersion <- function(response_matrix, covariate_matrix) {
  
  # perform precomputation and get dispersion one gene at a time
  response_indices <- seq_len(nrow(response_matrix))
  dispersions <- bplapply(response_indices, FUN = calculate_sceptre_dispersion_one_gene,
                          response_matrix = response_matrix,
                          covariate_matrix = covariate_matrix)
  
  # combine computed dispersions into one data frame
  dispersions <- rbindlist(dispersions)
  
  return(dispersions)
  
}

# run sceptre precomputation and extract dispersion for one gene
calculate_sceptre_dispersion_one_gene <- function(response_index, response_matrix,
                                                  covariate_matrix) {
  
  # report progress
  n_response <- paste0(" (", response_index, " of ", nrow(response_matrix), ")")
  response_id <- rownames(response_matrix)[[response_index]]
  message("  Calculating dispersion for response ", response_id, n_response)
  
  # perform Sceptre precomputation
  precomp <- sceptre:::perform_response_precomputation(
    expressions = response_matrix[response_index, ],
    covariate_matrix = covariate_matrix
  )
  
  # extract theta and create single row data frame as output
  dispersion <- data.table(response_id, expression_size = precomp$theta)
  
  return(dispersion)
  
}

# calculate gene expression baseline statistics (dispersion and mean expression) for all genes
calculate_expr_baseline <- function(sceptre_object, genes = NULL) {
  
  # process genes parameter or default to all genes if not specified
  if (!is.null(genes)) {
    genes <- intersect(genes, rownames(sceptre_object@response_matrix[[1]]))
    if (length(genes) == 0) {
      stop("None of the specified genes found in expression data", call. = FALSE)
    }
  } else {
    genes <- rownames(sceptre_object@response_matrix[[1]])
  }
  
  # run Sceptre precomputations to calculate dispersion for every gene
  dispersions <- calculate_sceptre_dispersion(
    response_matrix = sceptre_object@response_matrix[[1]][genes, ],
    covariate_matrix = sceptre_object@covariate_matrix
  )
  
  # calculate baseline mean expression for each gene
  expression_mean <- rowMeans(sceptre_object@response_matrix[[1]][genes, ])
  expression_mean <- data.table(response_id = names(expression_mean), expression_mean)
  
  # combine into one table with baseline gene expression statistics
  baseline_expression_stats <- merge(expression_mean, dispersions, by = "response_id")
  
  return(baseline_expression_stats)
  
}

## Perform power analysis --------------------------------------------------------------------------

# set up parallel backend if specified
if (opt$threads > 1) {
  message("Registering parallel backend with ", opt$threads, " threads")
  register(MulticoreParam(workers = opt$threads))
} else {
  message("Registering serial backend")
  register(SerialParam())
}

# load input MuData object
message("Loading input data...")
mudata <- readH5MU(opt$input)


## TODO: RUN SAME QC FILTERING AS INFERENCE ##


# get pairs to test table containing perturbation-response pairs to calculate power
pairs_to_test <- data.frame(metadata(mudata)$pairs_to_test, stringsAsFactors = FALSE)
discovery_pairs <- data.table(grna_target = pairs_to_test$intended_target_name,
                              response_id = pairs_to_test$gene_id)
discovery_pairs <- unique(discovery_pairs)  # make sure there are no duplicate pairs

# compute baseline expression stats (dispersion and mean expression) for all genes
message("Computing power analysis parameters:")
sceptre_object <- sceptreIGVF::convert_mudata_to_sceptre_object(mudata)
expression_stats <- calculate_expr_baseline(sceptre_object, genes = discovery_pairs$response_id)

# number of total cells in the experiment
num_total_cells <- ncol(mudata@ExperimentList$gene)

# get cells per guide RNA
cells_per_guide <- rowSums(assay(mudata@ExperimentList$guide, "guide_assignment"))
cells_per_guide <- data.table(grna_id = names(cells_per_guide), num_cells = cells_per_guide)

# get target for each guide
guide_metadata <- rowData(mudata@ExperimentList$guide)
guide_targets <- data.table(grna_id = rownames(guide_metadata),
                            grna_target = guide_metadata$intended_target_name)

# create cells_per_grna table for perturbplan
cells_per_grna <- merge(guide_targets, cells_per_guide, by = "grna_id")

# compute power for all discovery pairs using perturbplan
message("Computing power to detect ", opt$effect_size, " fold change effects...")
power_results <- compute_power_posthoc(
  num_total_cells = num_total_cells,
  cells_per_grna = cells_per_grna,
  baseline_expression_stats = expression_stats,
  discovery_pairs = discovery_pairs,
  control_group = opt$control_group,
  side = opt$side,
  n_nonzero_trt_thresh = opt$n_nonzero_trt_thresh,
  n_nonzero_cntrl_thresh = opt$n_nonzero_cntrl_thresh,
  fold_change_mean = opt$effect_size,
  fold_change_sd = opt$effect_size_sd,
  cutoff = opt$pvalue_cutoff
)

# add number of perturbed cells to individual power results
ind_power <- power_results$individual_power
n_cells <- cells_per_grna[, .(num_cells = sum(num_cells)), by = grna_target]
ind_power <- merge(ind_power, n_cells, by = "grna_target", all.x = TRUE)

# add gene expression stats (mean and dispersion) to individual power results
ind_power <- merge(ind_power, expression_stats, by = "response_id", all.x = TRUE)
setnames(ind_power, "expression_size", "expression_dispersion")

# add tested effect size to results list
power_results$effect_size <- opt$effect_size
power_results$individual_power <- ind_power

# add results to MuData object
message("Creating output MuData file...")
metadata(mudata)$power_results <- power_results

# write MuData object including power analysis results to new MuData file
writeH5MU(mudata, file = opt$output)


message("All done!")
