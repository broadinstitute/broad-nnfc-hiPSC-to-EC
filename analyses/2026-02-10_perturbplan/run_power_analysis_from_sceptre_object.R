#!/usr/bin/env Rscript

## Compute statistical power to detect a specific effect size on expression for all
## perturbation-gene pairs using a final sceptre object directly.

## Process input arguments
library(optparse)

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Path to input sceptre_object RDS file", metavar = "file"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Path to output RDS file with power analysis results", metavar = "file"),
  make_option("--expression_stats_output", type = "character", default = NULL,
              help = paste("Path to output TSV with per-gene fitted mean and theta.",
                           "Defaults to <output stem>_expression_stats.tsv"),
              metavar = "file"),
  make_option("--power_results_output", type = "character", default = NULL,
              help = paste("Path to output TSV with per-pair power results.",
                           "Defaults to <output stem>_individual_power.tsv"),
              metavar = "file"),
  make_option("--dropped_pairs_output", type = "character", default = NULL,
              help = paste("Path to output TSV with discovery pairs excluded before power analysis.",
                           "Defaults to <output stem>_dropped_pairs.tsv"),
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
                           "input sceptre discovery_result slot"),
              metavar = "real"),
  make_option("--control_group", type = "character", default = "complement",
              help = paste("String specifying the control group to use in the differential",
                           "expression analysis, either 'complement' or 'nt_cells'"),
              metavar = "character"),
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

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

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

default_dropped_pairs_output <- function(output) {
  paste0(tools::file_path_sans_ext(output), "_dropped_pairs.tsv")
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

if (is.null(opt$dropped_pairs_output)) {
  opt$dropped_pairs_output <- default_dropped_pairs_output(opt$output)
}

check_numeric_arg("effect_size", opt, lower = 0, lower_inclusive = FALSE)
check_numeric_arg("effect_size_sd", opt, lower = 0)
if (!is.null(opt$pvalue_cutoff)) {
  check_numeric_arg("pvalue_cutoff", opt, lower = 0, upper = 1,
                    lower_inclusive = FALSE, upper_inclusive = FALSE)
}

check_choice_arg("control_group", opt, choices = c("complement", "nt_cells"))
check_choice_arg("side", opt, choices = c("left", "right", "both"))
opt$n_nonzero_trt_thresh <- as.integer(check_integer_arg("n_nonzero_trt_thresh", opt, lower = 0))
opt$n_nonzero_cntrl_thresh <- as.integer(check_integer_arg("n_nonzero_cntrl_thresh", opt, lower = 0))

suppressPackageStartupMessages({
  library(Matrix)
  library(perturbplan)
  library(sceptre)
  library(data.table)
})

calculate_expression_stats <- function(sceptre_object, response_matrix, covariate_matrix,
                                       genes = NULL) {
  if (!is.null(genes)) {
    genes <- intersect(genes, rownames(response_matrix))
    if (length(genes) == 0) {
      stop("None of the specified genes found in expression data", call. = FALSE)
    }
  } else {
    genes <- rownames(response_matrix)
  }

  theta_vec <- vapply(
    sceptre_object@response_precomputations,
    function(x) x$theta,
    numeric(1)
  )
  theta_table <- data.table(response_id = names(theta_vec), expression_size = theta_vec)
  theta_table <- theta_table[response_id %in% genes]
  missing_theta_genes <- setdiff(genes, theta_table$response_id)
  if (length(missing_theta_genes) > 0) {
    stop(
      "Stored response_precomputations are missing theta for response_id values: ",
      paste(missing_theta_genes, collapse = ", "),
      call. = FALSE
    )
  }

  fitted_mean <- vapply(genes, function(gene_id) {
    precomp <- sceptre_object@response_precomputations[[gene_id]]
    if (is.null(precomp$fitted_coefs)) {
      stop("Stored response_precomputations are missing fitted_coefs for response_id ",
           gene_id, call. = FALSE)
    }
    mu <- exp(as.numeric(covariate_matrix %*% precomp$fitted_coefs))
    mean(mu)
  }, numeric(1))
  expression_stats <- data.table(
    response_id = genes,
    expression_mean = as.numeric(fitted_mean)
  )
  expression_stats <- merge(expression_stats, theta_table, by = "response_id")
  setnames(expression_stats, "expression_size", "expression_theta")
  expression_stats
}

get_package_versions <- function(packages) {
  data.table(
    package = packages,
    version = vapply(packages, function(package) {
      as.character(utils::packageVersion(package))
    }, character(1))
  )
}

get_cells_in_use <- function(sceptre_object) {
  cells_in_use <- sceptre_object@cells_in_use
  if (length(cells_in_use) == 0) {
    cells_in_use <- seq_len(ncol(sceptre_object@response_matrix[[1]]))
  }

  cells_in_use
}

derive_pvalue_cutoff_from_sceptre <- function(sceptre_object, discovery_pairs) {
  discovery_results <- as.data.table(sceptre_object@discovery_result)
  required_cols <- c("response_id", "grna_target", "p_value", "significant")
  if (!all(required_cols %in% colnames(discovery_results))) {
    stop("sceptre_object@discovery_result must contain columns: ",
         paste(required_cols, collapse = ", "), call. = FALSE)
  }

  discovery_pairs <- unique(as.data.table(discovery_pairs)[, .(grna_target, response_id)])
  discovery_results <- merge(
    discovery_results[, .(grna_target, response_id, p_value, significant)],
    discovery_pairs,
    by = c("grna_target", "response_id")
  )
  discovery_results <- discovery_results[!is.na(p_value)]
  if (nrow(discovery_results) == 0) {
    stop("No discovery rows with non-missing p_value found in sceptre_object@discovery_result",
         call. = FALSE)
  }

  significant_results <- discovery_results[significant %in% TRUE]
  if (nrow(significant_results) == 0) {
    stop("No significant discovery rows found in sceptre_object@discovery_result",
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
  output_dir <- dirname(output)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  fwrite(expression_stats, file = output, sep = "\t")
}

write_power_results <- function(ind_power, output) {
  output_dir <- dirname(output)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  fwrite(ind_power, file = output, sep = "\t")
}

write_dropped_pairs <- function(dropped_pairs, output) {
  output_dir <- dirname(output)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  fwrite(dropped_pairs, file = output, sep = "\t")
}

write_power_results_rds <- function(power_results, output) {
  output_dir <- dirname(output)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  saveRDS(power_results, file = output)
}

run_power_analysis <- function(baseline_expression_stats, discovery_pairs,
                               cells_per_grna, num_total_cells, opt) {
  message("Computing power using fitted expression mean...")
  power_results <- compute_power_posthoc(
    num_total_cells = num_total_cells,
    cells_per_grna = cells_per_grna,
    baseline_expression_stats = baseline_expression_stats,
    discovery_pairs = discovery_pairs,
    control_group = opt$control_group,
    side = opt$side,
    n_nonzero_trt_thresh = opt$n_nonzero_trt_thresh,
    n_nonzero_cntrl_thresh = opt$n_nonzero_cntrl_thresh,
    fold_change_mean = opt$effect_size,
    fold_change_sd = opt$effect_size_sd,
    cutoff = opt$pvalue_cutoff
  )

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

  ind_power <- merge(ind_power, baseline_expression_stats, by = "response_id", all.x = TRUE)
  if (nrow(ind_power) != num_ind_power_rows) {
    stop("Merging expression stats changed the number of individual power rows",
         call. = FALSE)
  }
  check_unique_pairs(ind_power, "individual power after merging expression stats")
  check_no_missing_values(ind_power, c("expression_mean", "expression_size"),
                          "individual power after merging expression stats")
  setnames(ind_power, "expression_size", "expression_theta")

  list(power_results = power_results, individual_power = ind_power)
}

extract_cells_per_grna <- function(sceptre_object, cells_in_use) {
  grna_assignment_matrix <- sceptre::get_grna_assignments(sceptre_object)
  if (is.null(grna_assignment_matrix)) {
    stop("sceptre::get_grna_assignments(sceptre_object) returned NULL", call. = FALSE)
  }

  grna_target_data_frame <- as.data.table(sceptre_object@grna_target_data_frame)
  required_cols <- c("grna_id", "grna_target")
  if (!all(required_cols %in% colnames(grna_target_data_frame))) {
    stop("sceptre_object@grna_target_data_frame must contain columns: ",
         paste(required_cols, collapse = ", "), call. = FALSE)
  }

  keep_guides <- !duplicated(grna_target_data_frame$grna_id)
  grna_target_data_frame <- grna_target_data_frame[keep_guides, .(grna_id, grna_target)]
  grna_assignment_matrix <- grna_assignment_matrix[keep_guides, cells_in_use, drop = FALSE]

  cells_per_guide <- data.table(
    grna_id = grna_target_data_frame$grna_id,
    num_cells = as.numeric(Matrix::rowSums(grna_assignment_matrix))
  )

  merge(grna_target_data_frame, cells_per_guide, by = "grna_id")
}

filter_discovery_pairs_to_available_targets <- function(discovery_pairs, available_grna_targets) {
  missing_discovery_targets <- setdiff(unique(discovery_pairs$grna_target), available_grna_targets)
  dropped_pairs <- discovery_pairs[0, .(grna_target, response_id)]
  if (length(missing_discovery_targets) > 0) {
    dropped_pairs <- discovery_pairs[
      grna_target %in% missing_discovery_targets,
      .(grna_target, response_id)
    ]
    dropped_pairs[, drop_reason := "missing_grna_target_in_assignments"]
    num_discarded_pairs <- nrow(dropped_pairs)
    message(
      "Dropping ", num_discarded_pairs,
      " discovery pairs whose grna_target is absent from the guide assignments (",
      length(missing_discovery_targets), " targets)."
    )
    discovery_pairs <- discovery_pairs[!(grna_target %in% missing_discovery_targets)]
  }

  if (nrow(discovery_pairs) == 0) {
    stop(
      "No discovery pairs remain after filtering to grna_target targets present in the guide assignments",
      call. = FALSE
    )
  }

  list(
    discovery_pairs = discovery_pairs,
    missing_discovery_targets = missing_discovery_targets,
    dropped_pairs = dropped_pairs
  )
}

filter_discovery_pairs_to_available_theta <- function(discovery_pairs, available_response_ids) {
  missing_response_ids <- setdiff(unique(discovery_pairs$response_id), available_response_ids)
  dropped_pairs <- discovery_pairs[0, .(grna_target, response_id)]
  if (length(missing_response_ids) > 0) {
    dropped_pairs <- discovery_pairs[
      response_id %in% missing_response_ids,
      .(grna_target, response_id)
    ]
    dropped_pairs[, drop_reason := "missing_theta_in_response_precomputations"]
    num_discarded_pairs <- nrow(dropped_pairs)
    message(
      "Dropping ", num_discarded_pairs,
      " discovery pairs whose response_id lacks stored theta in response_precomputations (",
      length(missing_response_ids), " genes)."
    )
    discovery_pairs <- discovery_pairs[!(response_id %in% missing_response_ids)]
  }

  if (nrow(discovery_pairs) == 0) {
    stop(
      "No discovery pairs remain after filtering to response_id values with stored theta in response_precomputations",
      call. = FALSE
    )
  }

  list(
    discovery_pairs = discovery_pairs,
    missing_response_ids = missing_response_ids,
    dropped_pairs = dropped_pairs
  )
}

pvalue_cutoff_source <- "command_line"
num_significant_sceptre_results <- NA_integer_
num_sceptre_results_for_cutoff <- NA_integer_

message("Loading input data...")
sceptre_object <- readRDS(opt$input)

cells_in_use <- get_cells_in_use(sceptre_object)
response_matrix <- sceptre_object@response_matrix[[1]][, cells_in_use, drop = FALSE]
covariate_matrix <- sceptre_object@covariate_matrix[cells_in_use, , drop = FALSE]

discovery_pairs <- unique(as.data.table(sceptre_object@discovery_pairs)[, .(grna_target, response_id)])
check_unique_pairs(discovery_pairs, "sceptre_object@discovery_pairs")

available_response_ids <- names(vapply(
  sceptre_object@response_precomputations,
  function(x) x$theta,
  numeric(1)
))
theta_filter <- filter_discovery_pairs_to_available_theta(
  discovery_pairs = discovery_pairs,
  available_response_ids = available_response_ids
)
discovery_pairs <- theta_filter$discovery_pairs
missing_theta_response_ids <- theta_filter$missing_response_ids
dropped_pairs_missing_theta <- theta_filter$dropped_pairs

if (is.null(opt$pvalue_cutoff)) {
  message("Deriving p-value cutoff from significant sceptre discovery results")
  pvalue_cutoff_info <- derive_pvalue_cutoff_from_sceptre(
    sceptre_object = sceptre_object,
    discovery_pairs = discovery_pairs
  )
  opt$pvalue_cutoff <- pvalue_cutoff_info$pvalue_cutoff
  pvalue_cutoff_source <- "sceptre_discovery_result_significant"
  num_significant_sceptre_results <- pvalue_cutoff_info$num_significant
  num_sceptre_results_for_cutoff <- pvalue_cutoff_info$num_tested
  message("Using p-value cutoff ", format(opt$pvalue_cutoff, digits = 17),
          " from the largest significant sceptre discovery p-value")
}

cells_per_grna <- extract_cells_per_grna(
  sceptre_object = sceptre_object,
  cells_in_use = cells_in_use
)
available_grna_targets <- unique(stats::na.omit(cells_per_grna$grna_target))
discovery_pair_filter <- filter_discovery_pairs_to_available_targets(
  discovery_pairs = discovery_pairs,
  available_grna_targets = available_grna_targets
)
discovery_pairs <- discovery_pair_filter$discovery_pairs
missing_discovery_targets <- discovery_pair_filter$missing_discovery_targets
dropped_pairs_missing_grna_targets <- discovery_pair_filter$dropped_pairs

dropped_pairs <- rbindlist(
  list(dropped_pairs_missing_theta, dropped_pairs_missing_grna_targets),
  use.names = TRUE,
  fill = TRUE
)
if (nrow(dropped_pairs) == 0) {
  dropped_pairs <- data.table(
    grna_target = character(),
    response_id = character(),
    drop_reason = character()
  )
}
message("Writing dropped-pairs log to ", opt$dropped_pairs_output)
write_dropped_pairs(dropped_pairs, opt$dropped_pairs_output)

message("Computing power analysis parameters:")
expression_stats <- calculate_expression_stats(
  sceptre_object = sceptre_object,
  response_matrix = response_matrix,
  covariate_matrix = covariate_matrix,
  genes = discovery_pairs$response_id
)
if (anyDuplicated(expression_stats$response_id)) {
  stop("Computed expression stats contain duplicate response_id rows", call. = FALSE)
}
message("Writing expression stats to ", opt$expression_stats_output)
write_expression_stats(expression_stats, opt$expression_stats_output)

num_total_cells <- ncol(response_matrix)

baseline_expression_stats <- copy(expression_stats[, .(response_id, expression_mean, expression_theta)])
setnames(baseline_expression_stats, "expression_theta", "expression_size")

power_analysis <- run_power_analysis(
  baseline_expression_stats = baseline_expression_stats,
  discovery_pairs = discovery_pairs,
  cells_per_grna = cells_per_grna,
  num_total_cells = num_total_cells,
  opt = opt
)

message("Writing power results to ", opt$power_results_output)
write_power_results(power_analysis$individual_power, opt$power_results_output)

power_results <- list(
  fitted_mean = power_analysis$power_results,
  expression_stats = expression_stats
)
power_results$effect_size <- opt$effect_size
power_results$run_parameters <- data.table(
  run_timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  input = normalizePath(opt$input, mustWork = FALSE),
  output = normalizePath(opt$output, mustWork = FALSE),
  expression_stats_output = normalizePath(opt$expression_stats_output, mustWork = FALSE),
  power_results_output = normalizePath(opt$power_results_output, mustWork = FALSE),
  dropped_pairs_output = normalizePath(opt$dropped_pairs_output, mustWork = FALSE),
  command_args = paste(commandArgs(trailingOnly = TRUE), collapse = " "),
  effect_size = opt$effect_size,
  effect_size_sd = opt$effect_size_sd,
  pvalue_cutoff = opt$pvalue_cutoff,
  pvalue_cutoff_source = pvalue_cutoff_source,
  num_significant_sceptre_results = num_significant_sceptre_results,
  num_sceptre_results_for_cutoff = num_sceptre_results_for_cutoff,
  control_group = opt$control_group,
  side = opt$side,
  n_nonzero_trt_thresh = opt$n_nonzero_trt_thresh,
  n_nonzero_cntrl_thresh = opt$n_nonzero_cntrl_thresh,
  expression_stats_source = "fitted_from_input_sceptre_object",
  cells_in_use_source = if (length(sceptre_object@cells_in_use) == 0) "all_cells" else "sceptre_object@cells_in_use",
  num_total_cells = num_total_cells,
  num_cells_in_use = length(cells_in_use),
  num_discovery_pairs = nrow(discovery_pairs),
  num_response_genes = length(unique(discovery_pairs$response_id)),
  num_grna_targets = length(unique(discovery_pairs$grna_target)),
  num_dropped_pairs = nrow(dropped_pairs),
  num_dropped_pairs_missing_grna_target = nrow(dropped_pairs_missing_grna_targets),
  num_dropped_pairs_missing_theta = nrow(dropped_pairs_missing_theta),
  num_missing_grna_targets = length(missing_discovery_targets),
  num_missing_theta_response_ids = length(missing_theta_response_ids),
  low_moi = sceptre_object@low_moi,
  r_version = R.version.string
)
power_results$package_versions <- get_package_versions(
  c("Matrix", "perturbplan", "sceptre", "data.table")
)
power_results$missing_grna_targets <- missing_discovery_targets
power_results$missing_theta_response_ids <- missing_theta_response_ids
power_results$dropped_pairs <- dropped_pairs

message("Writing power results RDS to ", opt$output)
write_power_results_rds(power_results, opt$output)

message("All done!")