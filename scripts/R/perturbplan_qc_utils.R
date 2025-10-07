#!/usr/bin/env Rscript
# Author: Eugenio Mattei (Broad Institute of MIT and Harvard)

# ------------------------------------------------------------------------
# perturbplan Utilities for QC
# ------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(BiocParallel)
  library(sceptre)
  library(data.table)
})



#' Calculate SCEPTRE dispersion for all genes
#'
#' @param response_matrix Matrix of gene expression (genes x cells)
#' @param covariate_matrix Matrix of covariates (cells x covariates)
#' @return data.table with response_id and expression_size for all genes
calculate_sceptre_dispersion <- function(response_matrix, covariate_matrix) {
  
  # Compute dispersion one gene at a time using bplapply
  response_indices <- seq_len(nrow(response_matrix))
  dispersions <- BiocParallel::bplapply(
    response_indices,
    FUN = .calculate_sceptre_dispersion_helper,
    response_matrix = response_matrix,
    covariate_matrix = covariate_matrix
    )
  
  # combine computed dispersions into one data frame
  dispersions <- rbindlist(dispersions)
  
  return(dispersions)
  
}

#' Helper function: Calculate SCEPTRE dispersion for a single gene
#'
#' This is a helper function to be used with bplapply to calculate SCEPTRE
#' dispersion for one gene at a time.
#'
#' @param response_index Integer index of the gene (row) in response_matrix
#' @param response_matrix Matrix of gene expression (genes x cells)
#' @param covariate_matrix Matrix of covariates (cells x covariates)
#' @return data.table with response_id and expression_size (dispersion)
.calculate_sceptre_dispersion_helper <- function(response_index, response_matrix,
                                                  covariate_matrix) {
  
  # report progress
  n_response <- sprintf(" (%d of %d)", response_index, nrow(response_matrix))
  response_id <- rownames(response_matrix)[[response_index]]
  message("  Calculating dispersion for response ", response_id, n_response)
  
  # Use a Poisson GLM to get 𝜇 for the gene quickly, and then profile for θ.
  # (glm.nb does internally but with fewer iterations). Downstream, 𝜇 and 𝜃
  # are used to compute the mean/SD of a Z-like test statistic for power calculations
  # Large θ → Poisson-like (little overdispersion).
  # Small θ → strong overdispersion
  # Check compute_distribution_teststat(...)
  # “dispersion” 𝛼 it’s 𝛼 = 1/𝜃.
  # Returns:
  # fitted_coefs: the Poisson 𝛽 (β not refitted under NB; fast approximation)
  # theta: the NB size/overdispersion used later to set test-stat variability
  precomp <- sceptre:::perform_response_precomputation(
    expressions = response_matrix[response_index, ],
    covariate_matrix = covariate_matrix
  )
  
  # extract theta and create single row data frame as output
  dispersion <- data.table(response_id, expression_size = precomp$theta)
  
  return(dispersion)
  
}

#' Calculate baseline gene expression statistics (mean and dispersion)
#'
#' @param sceptre_object SCEPTRE object containing response_matrix and covariate_matrix
#' @param genes Optional vector of gene names to restrict analysis
#' @return data.table with response_id, expression_mean, and expression_size
calculate_expr_baseline <- function(sceptre_object, genes = NULL) {
  
  # Process genes parameter or default to all genes if not specified
  all_genes <- rownames(sceptre_object@response_matrix[[1]])
  genes <- if (!is.null(genes)) intersect(genes, all_genes) else all_genes
  if (length(genes) == 0) {
    stop("None of the specified genes found in expression data", call. = FALSE)
  }
  
  # Calculate dispersion for every gene.
  dispersions <- calculate_sceptre_dispersion(
    response_matrix = sceptre_object@response_matrix[[1]][genes, ],
    covariate_matrix = sceptre_object@covariate_matrix
  )
  
  # calculate baseline mean expression for each gene
  expression_mean <- Matrix::rowMeans(sceptre_object@response_matrix[[1]][genes, ])
  expression_mean <- data.table(response_id = names(expression_mean), expression_mean)
  
  # combine into one table with baseline gene expression statistics
  baseline_expression_stats <- merge(expression_mean, dispersions, by = "response_id")
  
  return(baseline_expression_stats)
  
}

get_avg_num_cells_per_grna <- function(sceptre_object){
  return(mean(Matrix::rowSums(methods::as(sceptre::get_grna_assignments(sceptre_object, apply_cellwise_qc = FALSE),"dsparseMatrix"))))
}
