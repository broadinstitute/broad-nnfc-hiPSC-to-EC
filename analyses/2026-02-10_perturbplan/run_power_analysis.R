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
  make_option("--expression_stats_output", type = "character", default = NULL,
              help = paste("Path to output TSV with computed per-gene expression mean and theta.",
                           "Defaults to <output stem>_expression_stats.tsv"),
              metavar = "file"),
  make_option("--power_results_output", type = "character", default = NULL,
              help = paste("Path to output TSV with per-pair power results.",
                           "Defaults to <output stem>_individual_power.tsv"),
              metavar = "file"),
  make_option(c("-e", "--effect_size"), type = "double", default = NULL,
              help = "Fold change effect size on expression, for which power will be computed",
              metavar = "real"),
  make_option(c("-s", "--effect_size_sd"), type = "double", default = 0.13,
              help = "Effect size standard deviation to model guide-guide variability",
              metavar = "real"),
  make_option(c("-p", "--pvalue_cutoff"), type = "double", default = NULL,
              help = paste("Nominal p-value cutoff used in power analysis.",
                           "If omitted, this is set to the largest discovery",
                           "p-value among rows marked significant in the",
                           "input MuData test_results metadata"),
              metavar = "real"),
  make_option(c("-t", "--threads"), type = "double", default = 1,
              help = "Number of CPU threads to use for parallelization", metavar = "integer"),
  make_option("--control_group", type = "character", default = "complement",
              help = paste("String specifying the control group to use in the differential",
                "expression analysis, either 'complement' or 'nt_cells'"), metavar = "character"),
  make_option("--side", type = "character", default = "both",
              help = "The sidedness of the test, one of 'left', 'right', or 'both'",
              metavar = "character"),
  make_option("--n_nonzero_trt_thresh", type = "double", default = 7,
              help = "Number of nonzero treatment cells a pair must contain for it to be retained",
              metavar = "integer"),
  make_option("--n_nonzero_cntrl_thresh", type = "double", default = 7,
              help = "Number of nonzero control cells a pair must contain for it to be retained",
              metavar = "integer")

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

stop_bad_arg <- function(arg, message) {
  stop("--", arg, " ", message, call. = FALSE)
}

default_expression_stats_output <- function(output) {
  paste0(tools::file_path_sans_ext(output), "_expression_stats.tsv")
}

default_power_results_output <- function(output) {
  paste0(tools::file_path_sans_ext(output), "_individual_power.tsv")
}

check_numeric_arg <- function(arg, opt, lower = -Inf, upper = Inf,
                              lower_inclusive = TRUE, upper_inclusive = TRUE) {
  value <- opt[[arg]]
  if (!is.numeric(value) || length(value) != 1 || is.na(value) || !is.finite(value)) {
    stop_bad_arg(arg, "must be a finite number")
  }

  lower_failed <- if (lower_inclusive) value < lower else value <= lower
  if (lower_failed) {
    operator <- if (lower_inclusive) ">=" else ">"
    stop_bad_arg(arg, paste("must be", operator, lower))
  }

  upper_failed <- if (upper_inclusive) value > upper else value >= upper
  if (upper_failed) {
    operator <- if (upper_inclusive) "<=" else "<"
    stop_bad_arg(arg, paste("must be", operator, upper))
  }

  invisible(value)
}

check_integer_arg <- function(arg, opt, lower = -Inf) {
  value <- check_numeric_arg(arg, opt, lower = lower)
  if (value %% 1 != 0) {
    stop_bad_arg(arg, "must be an integer")
  }

  invisible(value)
}

check_choice_arg <- function(arg, opt, choices) {
  value <- opt[[arg]]
  if (!is.character(value) || length(value) != 1 || is.na(value) ||
      !(value %in% choices)) {
    stop_bad_arg(arg, paste("must be one of:", paste(choices, collapse = ", ")))
  }

  invisible(value)
}

# check that all required parameters are provided
required_args <- c("input", "output", "effect_size")
for (arg in required_args) {
  check_required_args(arg, opt = opt, opt_parser = opt_parser)
}

if (is.null(opt$expression_stats_output)) {
  opt$expression_stats_output <- default_expression_stats_output(opt$output)
}

if (is.null(opt$power_results_output)) {
  opt$power_results_output <- default_power_results_output(opt$output)
}

# check that argument values are valid before loading data or running expensive steps
check_numeric_arg("effect_size", opt, lower = 0, lower_inclusive = FALSE)
check_numeric_arg("effect_size_sd", opt, lower = 0)
if (!is.null(opt$pvalue_cutoff)) {
  check_numeric_arg("pvalue_cutoff", opt, lower = 0, upper = 1,
                    lower_inclusive = FALSE, upper_inclusive = FALSE)
}

opt$threads <- as.integer(check_integer_arg("threads", opt, lower = 1))
check_choice_arg("control_group", opt, choices = c("complement", "nt_cells"))
check_choice_arg("side", opt, choices = c("left", "right", "both"))
opt$n_nonzero_trt_thresh <- as.integer(check_integer_arg("n_nonzero_trt_thresh", opt, lower = 0))
opt$n_nonzero_cntrl_thresh <- as.integer(check_integer_arg("n_nonzero_cntrl_thresh", opt, lower = 0))

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

# run sceptre precomputation on all genes to get negative binomial theta values
calculate_sceptre_theta <- function(response_matrix, covariate_matrix) {

  # perform precomputation and get theta one gene at a time
  response_indices <- seq_len(nrow(response_matrix))
  thetas <- bplapply(response_indices, FUN = calculate_sceptre_theta_one_gene,
                     response_matrix = response_matrix,
                     covariate_matrix = covariate_matrix)

  # combine computed theta values into one data frame
  thetas <- rbindlist(thetas)

  return(thetas)

}

perform_sceptre_response_precomputation <- function(expressions, covariate_matrix,
                                                    response_id) {
  precompute_response <- getFromNamespace("perform_response_precomputation", "sceptre")
  tryCatch(
    precompute_response(
      expressions = expressions,
      covariate_matrix = covariate_matrix
    ),
    error = function(e) {
      stop("Failed to calculate Sceptre theta for response ", response_id, ": ",
           conditionMessage(e), call. = FALSE)
    }
  )
}

# run sceptre precomputation and extract theta for one gene
calculate_sceptre_theta_one_gene <- function(response_index, response_matrix,
                                             covariate_matrix) {

  # report progress
  n_response <- paste0(" (", response_index, " of ", nrow(response_matrix), ")")
  response_id <- rownames(response_matrix)[[response_index]]
  message("  Calculating theta for response ", response_id, n_response)

  # perform Sceptre precomputation
  precomp <- perform_sceptre_response_precomputation(
    expressions = as.numeric(response_matrix[response_index, , drop = TRUE]),
    covariate_matrix = covariate_matrix,
    response_id = response_id
  )

  # extract theta as the negative binomial size expected by perturbplan
  theta <- data.table(response_id, expression_size = precomp$theta)

  return(theta)

}

# calculate gene expression baseline statistics (theta and mean expression) for all genes
calculate_expr_baseline <- function(sceptre_object, genes = NULL) {

  response_matrix <- sceptre_object@response_matrix[[1]]

  # process genes parameter or default to all genes if not specified
  if (!is.null(genes)) {
    genes <- intersect(genes, rownames(response_matrix))
    if (length(genes) == 0) {
      stop("None of the specified genes found in expression data", call. = FALSE)
    }
  } else {
    genes <- rownames(response_matrix)
  }

  response_matrix <- response_matrix[genes, , drop = FALSE]

  # run Sceptre precomputations to calculate theta for every gene
  thetas <- calculate_sceptre_theta(
    response_matrix = response_matrix,
    covariate_matrix = sceptre_object@covariate_matrix
  )

  # calculate baseline mean expression for each gene
  expression_mean <- rowMeans(response_matrix)
  expression_mean <- data.table(response_id = names(expression_mean), expression_mean)

  # combine into one table with baseline gene expression statistics
  baseline_expression_stats <- merge(expression_mean, thetas, by = "response_id")

  return(baseline_expression_stats)

}

get_package_versions <- function(packages) {
  data.table(
    package = packages,
    version = vapply(packages, function(package) {
      as.character(utils::packageVersion(package))
    }, character(1))
  )
}

derive_pvalue_cutoff_from_mudata <- function(mudata, discovery_pairs) {
  test_results <- as.data.table(
    data.frame(metadata(mudata)$test_results, stringsAsFactors = FALSE)
  )
  required_cols <- c("gene_id", "intended_target_name", "p_value", "pair_type",
                     "significant")
  if (!all(required_cols %in% colnames(test_results))) {
    stop("metadata(mudata)$test_results must contain columns: ",
         paste(required_cols, collapse = ", "), call. = FALSE)
  }

  discovery_pairs <- unique(as.data.table(discovery_pairs)[, .(grna_target, response_id)])
  discovery_results <- test_results[pair_type == "discovery", .(
    grna_target = intended_target_name,
    response_id = gene_id,
    p_value,
    significant
  )]
  discovery_results <- merge(
    discovery_results,
    discovery_pairs,
    by = c("grna_target", "response_id")
  )
  discovery_results <- discovery_results[!is.na(p_value)]
  if (nrow(discovery_results) == 0) {
    stop("No discovery rows with non-missing p_value found in metadata(mudata)$test_results",
         call. = FALSE)
  }

  significant_results <- discovery_results[significant %in% TRUE]
  if (nrow(significant_results) == 0) {
    stop("No significant discovery rows found in metadata(mudata)$test_results",
         call. = FALSE)
  }

  pvalue_cutoff <- max(significant_results$p_value)
  if (!is.finite(pvalue_cutoff) || pvalue_cutoff <= 0 || pvalue_cutoff >= 1) {
    stop("Derived p-value cutoff must be strictly between 0 and 1, got: ",
         pvalue_cutoff, call. = FALSE)
  }

  list(
    pvalue_cutoff = pvalue_cutoff,
    num_significant = nrow(significant_results),
    num_tested = nrow(discovery_results)
  )
}

check_unique_pairs <- function(data, context) {
  pair_keys <- paste(data$grna_target, data$response_id, sep = "\r")
  duplicate_keys <- unique(pair_keys[duplicated(pair_keys)])
  if (length(duplicate_keys) > 0) {
    stop(context, " contains duplicate grna_target/response_id rows", call. = FALSE)
  }

  invisible(NULL)
}

check_no_missing_values <- function(data, columns, context) {
  data <- as.data.frame(data)
  missing_columns <- columns[vapply(data[, columns, drop = FALSE], anyNA, logical(1))]
  if (length(missing_columns) > 0) {
    stop(context, " has missing values in columns: ",
         paste(missing_columns, collapse = ", "), call. = FALSE)
  }

  invisible(NULL)
}

write_expression_stats <- function(expression_stats, output) {
  expression_stats_output <- copy(expression_stats)
  setnames(expression_stats_output, "expression_size", "expression_theta")

  output_dir <- dirname(output)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  fwrite(expression_stats_output, file = output, sep = "\t")
}

write_power_results <- function(ind_power, output) {
  output_dir <- dirname(output)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  fwrite(ind_power, file = output, sep = "\t")
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

pvalue_cutoff_source <- "command_line"
num_significant_mudata_results <- NA_integer_
num_mudata_results_for_cutoff <- NA_integer_

# load input MuData object
message("Loading input data...")
mudata <- readH5MU(opt$input)
mudata_experiments <- MultiAssayExperiment::experiments(mudata)
required_experiments <- c("gene", "guide")
if (!all(required_experiments %in% names(mudata_experiments))) {
  stop("Input MuData must contain experiments: ",
       paste(required_experiments, collapse = ", "), call. = FALSE)
}

gene_experiment <- mudata_experiments[["gene"]]
guide_experiment <- mudata_experiments[["guide"]]


## TODO: RUN SAME QC FILTERING AS INFERENCE ##


# get discovery perturbation-response pairs to calculate power
pairs_to_test <- data.frame(metadata(mudata)$pairs_to_test, stringsAsFactors = FALSE)
required_pair_cols <- c("intended_target_name", "gene_id", "pair_type")
if (!all(required_pair_cols %in% colnames(pairs_to_test))) {
  stop("metadata(mudata)$pairs_to_test must contain columns: ",
       paste(required_pair_cols, collapse = ", "), call. = FALSE)
}

pairs_to_test <- pairs_to_test[pairs_to_test$pair_type == "discovery", ]
if (nrow(pairs_to_test) == 0) {
  stop("No discovery pairs found in metadata(mudata)$pairs_to_test", call. = FALSE)
}

discovery_pairs <- data.table(grna_target = pairs_to_test$intended_target_name,
                              response_id = pairs_to_test$gene_id)
discovery_pairs <- unique(discovery_pairs)  # make sure there are no duplicate pairs

if (is.null(opt$pvalue_cutoff)) {
  message("Deriving p-value cutoff from significant MuData discovery test results")
  pvalue_cutoff_info <- derive_pvalue_cutoff_from_mudata(
    mudata = mudata,
    discovery_pairs = discovery_pairs
  )
  opt$pvalue_cutoff <- pvalue_cutoff_info$pvalue_cutoff
  pvalue_cutoff_source <- "mudata_test_results_significant"
  num_significant_mudata_results <- pvalue_cutoff_info$num_significant
  num_mudata_results_for_cutoff <- pvalue_cutoff_info$num_tested
  message("Using p-value cutoff ", format(opt$pvalue_cutoff, digits = 17),
          " from the largest significant MuData discovery p-value")
}

# compute baseline expression stats (theta and mean expression) for all genes
message("Computing power analysis parameters:")
sceptre_object <- sceptreIGVF::convert_mudata_to_sceptre_object(mudata)
expression_stats <- calculate_expr_baseline(sceptre_object, genes = discovery_pairs$response_id)
if (anyDuplicated(expression_stats$response_id)) {
  stop("Computed expression stats contain duplicate response_id rows", call. = FALSE)
}
message("Writing expression stats to ", opt$expression_stats_output)
write_expression_stats(expression_stats, opt$expression_stats_output)

# number of total cells in the experiment
num_total_cells <- ncol(gene_experiment)

# get cells per guide RNA
cells_per_guide <- rowSums(SummarizedExperiment::assay(guide_experiment, "guide_assignment"))
cells_per_guide <- data.table(grna_id = names(cells_per_guide), num_cells = cells_per_guide)

# get target for each guide
guide_metadata <- SummarizedExperiment::rowData(guide_experiment)
guide_targets <- data.table(grna_id = rownames(guide_metadata),
                            grna_target = guide_metadata$intended_target_name)

# create cells_per_grna table for perturbplan
cells_per_grna <- merge(guide_targets, cells_per_guide, by = "grna_id")

available_grna_targets <- unique(stats::na.omit(cells_per_grna$grna_target))
missing_discovery_targets <- setdiff(unique(discovery_pairs$grna_target), available_grna_targets)
if (length(missing_discovery_targets) > 0) {
  num_discarded_pairs <- nrow(discovery_pairs[grna_target %in% missing_discovery_targets])
  message(
    "Dropping ", num_discarded_pairs,
    " discovery pairs whose grna_target is absent from the guide experiment (",
    length(missing_discovery_targets), " targets)."
  )
  discovery_pairs <- discovery_pairs[!(grna_target %in% missing_discovery_targets)]
}
if (nrow(discovery_pairs) == 0) {
  stop(
    "No discovery pairs remain after filtering to grna_target targets present in the guide experiment",
    call. = FALSE
  )
}

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
check_unique_pairs(ind_power, "power_results$individual_power")
num_ind_power_rows <- nrow(ind_power)

n_cells <- cells_per_grna[, .(num_cells = sum(num_cells)), by = grna_target]
ind_power <- merge(ind_power, n_cells, by = "grna_target", all.x = TRUE)
if (nrow(ind_power) != num_ind_power_rows) {
  stop("Merging guide-cell counts changed the number of individual power rows",
       call. = FALSE)
}
check_unique_pairs(ind_power, "individual power after merging guide-cell counts")
check_no_missing_values(ind_power, "num_cells",
                        "individual power after merging guide-cell counts")

# add gene expression stats (mean and theta) to individual power results
ind_power <- merge(ind_power, expression_stats, by = "response_id", all.x = TRUE)
if (nrow(ind_power) != num_ind_power_rows) {
  stop("Merging expression stats changed the number of individual power rows",
       call. = FALSE)
}
check_unique_pairs(ind_power, "individual power after merging expression stats")
check_no_missing_values(ind_power, c("expression_mean", "expression_size"),
                        "individual power after merging expression stats")
setnames(ind_power, "expression_size", "expression_theta")
message("Writing individual power results to ", opt$power_results_output)
write_power_results(ind_power, opt$power_results_output)

# add run provenance to results list
power_results$effect_size <- opt$effect_size
power_results$run_parameters <- data.table(
  run_timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  input = normalizePath(opt$input, mustWork = FALSE),
  output = normalizePath(opt$output, mustWork = FALSE),
  expression_stats_output = normalizePath(opt$expression_stats_output, mustWork = FALSE),
  power_results_output = normalizePath(opt$power_results_output, mustWork = FALSE),
  command_args = paste(commandArgs(trailingOnly = TRUE), collapse = " "),
  effect_size = opt$effect_size,
  effect_size_sd = opt$effect_size_sd,
  pvalue_cutoff = opt$pvalue_cutoff,
  pvalue_cutoff_source = pvalue_cutoff_source,
  num_significant_mudata_results = num_significant_mudata_results,
  num_mudata_results_for_cutoff = num_mudata_results_for_cutoff,
  control_group = opt$control_group,
  side = opt$side,
  n_nonzero_trt_thresh = opt$n_nonzero_trt_thresh,
  n_nonzero_cntrl_thresh = opt$n_nonzero_cntrl_thresh,
  threads = opt$threads,
  expression_stats_source = "computed_from_input_mudata",
  num_total_cells = num_total_cells,
  num_discovery_pairs = nrow(discovery_pairs),
  num_response_genes = length(unique(discovery_pairs$response_id)),
  num_grna_targets = length(unique(discovery_pairs$grna_target)),
  r_version = R.version.string
)
power_results$package_versions <- get_package_versions(
  c("MuData", "MultiAssayExperiment", "BiocParallel", "perturbplan",
    "sceptre", "sceptreIGVF", "data.table", "SummarizedExperiment")
)
power_results$individual_power <- ind_power

# add results to MuData object
message("Creating output MuData file...")
metadata(mudata)$power_results <- power_results

# write MuData object including power analysis results to new MuData file
writeH5MU(mudata, file = opt$output)


message("All done!")
