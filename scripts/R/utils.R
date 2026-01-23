# Check if ggtext is installed
library(ggtext)

# ---- Color Palettes ----

#' Qualitative color palette (4 colors)
#'
#' Useful for categorical data visualization.
palette_4_qual <- c("#D95F02", "#1B9E77", "#7570B3", "#E6AB02")

palette_4_qual_alt <- c("#D95F02", "#1B9E77", "#7570B3", "#E76F51")

#' Softer qualitative color palette (4 colors)
#'
#' Alternative palette with softer tones.
palette_4_qual_soft <- c("#264653", "#2A9D8F", "#E9C46A", "#E76F51")

# ---- Functions ----

#' Read 10x Genomics HDF5 file with RNA and optional CRISPR/gRNA modalities
#'
#' Reads a 10x h5 file, handling multiple modalities and nested genome lists.
#' Returns a Seurat object with RNA counts and, if present, a gRNA assay.
#'
#' @param h5 Path to the 10x HDF5 file.
#' @param sample_name Name for the Seurat project/sample.
#' @return Seurat object with RNA and optional gRNA assay.
#' @examples
#' obj <- read_10x_rna_grna("sample.h5", "Sample1")
read_10x_rna_grna <- function(h5, sample_name, load_gRNA=TRUE) {
  # Read the HDF5 file using Seurat's Read10X_h5
  x <- Seurat::Read10X_h5(h5, use.names = TRUE)
  
  # If x is nested by genome (list of lists), flatten one level
  if (is.list(x) && !inherits(x, "dgCMatrix") && all(vapply(x, is.list, logical(1)))) {
    x <- unlist(x, recursive = FALSE)
  }
  
  # Handle single-modality file (just RNA)
  if (inherits(x, "dgCMatrix")) {
    obj <- Seurat::CreateSeuratObject(counts = x, project = sample_name)
    return(obj)
  }
  
  # Identify modality names case-insensitively
  nms <- names(x)
  
  # Helper function to pick the correct modality name
  pick_name <- function(cands) {
    hit <- nms[match(tolower(cands), tolower(nms))]
    hit[!is.na(hit)][1]
  }
  
  # Try to find RNA and gRNA modalities
  rna_name  <- pick_name(c("Gene Expression","gene expression","GeneExpression","RNA"))
  grna_name <- pick_name(c("CRISPR Guide Capture","CRISPR Guides","CRISPR","Guide Capture"))
  
  if (is.na(rna_name)) stop("Could not find an RNA modality in: ", paste(nms, collapse=", "))
  
  # Create Seurat object with RNA counts
  obj <- Seurat::CreateSeuratObject(counts = x[[rna_name]], project = sample_name)
  
  # If gRNA modality is present, add as assay
  if (!is.na(grna_name) && load_gRNA) {
    obj[["gRNA"]] <- Seurat::CreateAssayObject(counts = x[[grna_name]])
  } else {
    message("No CRISPR/gRNA modality found in: ", h5)
  }
  obj <- Seurat::RenameCells(obj, add.cell.id = sample_name)
  
  obj
}

read_10x_grna <- function(h5, sample_name) {
  # Read the HDF5 file using Seurat's Read10X_h5
  x <- Seurat::Read10X_h5(h5, use.names = TRUE)
  
  # If x is nested by genome (list of lists), flatten one level
  if (is.list(x) && !inherits(x, "dgCMatrix") && all(vapply(x, is.list, logical(1)))) {
    x <- unlist(x, recursive = FALSE)
  }
  
  # Identify modality names case-insensitively
  nms <- names(x)
  
  # Helper function to pick the correct modality name
  pick_name <- function(cands) {
    hit <- nms[match(tolower(cands), tolower(nms))]
    hit[!is.na(hit)][1]
  }
  
  # Try to find gRNA modality
  grna_name <- pick_name(c("CRISPR Guide Capture","CRISPR Guides","CRISPR","Guide Capture"))
  
  if (is.na(grna_name)) stop("Could not find a gRNA modality in: ", paste(nms, collapse=", "))
  
  obj <- Seurat::CreateSeuratObject(counts = x[[grna_name]], project = sample_name, assay = "gRNA")
  
  obj <- Seurat::RenameCells(obj, add.cell.id = sample_name)
  
  obj
}



# Create function to save ggplot
save_ggplot <- function(plot, filename, width = 8, height = 6, dpi = 300) {
  ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi
  )
}


# ---- plotting + saving helpers ----
# Save pdf and png
save_both <- function(plot, out_dir, stem, width = 8, height = 6, dpi = 300) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(out_dir, paste0(stem, ".png")), plot = plot, width = width, height = height, dpi = dpi)
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), plot = plot, width = width, height = height, dpi = dpi)
  invisible(plot)
}

# Plot histogram for a value but mark the median
plot_hist_with_group_medians <- function(df, value_col, group_col = "day_factor",
                                         title, xlab, bins = 50) {
  stopifnot(value_col %in% names(df), group_col %in% names(df))
  
  # histogram data (global, used to place labels consistently)
  hist_data <- ggplot_build(
    ggplot(df, aes(x = .data[[value_col]])) +
      geom_histogram(bins = bins, position = "identity")
  )$data[[1]]
  y_label <- max(hist_data$count, na.rm = TRUE) * 0.5
  
  median_stats <- dplyr::summarise(
    dplyr::group_by(df, .data[[group_col]]),
    median_val = median(.data[[value_col]], na.rm = TRUE)
  )
  
  ggplot(df, aes(x = .data[[value_col]])) +
    geom_histogram(bins = bins, color = "black", position = "identity", alpha = 0.6) +
    geom_vline(
      data = median_stats,
      aes(xintercept = median_val),
      linetype = "dashed",
      color = "red",
      linewidth = 1
    ) +
    geom_text(
      data = median_stats,
      aes(x = median_val, y = y_label, label = paste0("Median: ", round(median_val, 2))),
      hjust = -0.5,
      vjust = 0.5,
      show.legend = FALSE,
      color = "black"
    ) +
    labs(title = title, x = xlab, y = "Barcode Count") +
    facet_grid(rows = vars(.data[[group_col]])) +
    theme_minimal(base_size = 14)
}

make_threshold_labels <- function(breaks, threshold, color, suffix = "") {
  labels <- as.character(breaks)
  labels[breaks == threshold] <- sprintf("<span style='color:%s;'>%s%s</span>", color, threshold, suffix)
  labels
}

# This replaces the Seurat::FeatureScatter function
plot_scatter_thresholds <- function(df,
                                    x = "nCount_RNA",
                                    y,
                                    title,
                                    subtitle = NULL,
                                    x_thresh = 1000,
                                    y_thresh = 50,
                                    x_max = 30000,
                                    x_breaks = seq(0, 30000, by = 5000),
                                    y_breaks = c(0, 20, 40, 50, 60, 80, 100),
                                    # --- marginal options ---
                                    add_marginals = TRUE,
                                    marginal_type = c("histogram", "density"),
                                    marginal_bins = NULL,
                                    min_bins = 30,
                                    max_bins = 100,
                                    marginal_fill = "grey80",
                                    marginal_alpha = 0.8) {
  
  stopifnot(x %in% names(df), y %in% names(df))
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
  if (!requireNamespace("ggtext", quietly = TRUE)) {
    stop("Package 'ggtext' is required. Install it with: install.packages('ggtext')")
  }
  
  if (isTRUE(add_marginals) && !requireNamespace("ggExtra", quietly = TRUE)) {
    stop("Package 'ggExtra' is required for marginals. Install it with: install.packages('ggExtra')")
  }
  
  marginal_type <- match.arg(marginal_type)
  
  auto_bins <- function(x, min_bins = 30, max_bins = 100) {
    n <- sum(!is.na(x))
    bins <- floor(n^(1/3))
    max(min_bins, min(max_bins, bins))
  }
  
  # Choose bins automatically unless user provided them
  bins_use <- if (is.null(marginal_bins)) {
    max(
      auto_bins(df[[x]], min_bins, max_bins),
      auto_bins(df[[y]], min_bins, max_bins)
    )
  } else {
    marginal_bins
  }
  
  x_breaks <- sort(unique(c(x_breaks, x_thresh)))
  y_breaks <- sort(unique(c(y_breaks, y_thresh)))
  
  x_labels <- as.character(x_breaks)
  y_labels <- as.character(y_breaks)
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x]], y = .data[[y]])) +
    ggplot2::geom_point(size = 0.3, alpha = 0.5) +
    ggplot2::geom_hline(yintercept = y_thresh, linetype = "dashed", color = "red") +
    ggplot2::geom_vline(xintercept = x_thresh, linetype = "dashed", color = "blue") +
    ggplot2::scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels,
      limits = c(0, x_max)
    ) +
    ggplot2::scale_y_continuous(
      breaks = y_breaks,
      labels = y_labels,
      limits = c(0, 100)
    ) +
    ggplot2::labs(
      title = title,
      subtitle = paste0(subtitle," bin_size = ", bins_use),
      x = x,
      y = y
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.text.x = ggtext::element_markdown(),
      axis.text.y = ggtext::element_markdown()
    )
  
  if (!isTRUE(add_marginals)) return(p)
  
  ggExtra::ggMarginal(
    p,
    type = marginal_type,
    bins = bins_use,
    fill = marginal_fill,
    alpha = marginal_alpha,
    colour = NA
  )
}

  # ---- target UMIs computation helper (Seurat v5 layers) ----

compute_percent_target_umis <- function(seurat_obj, target_genes, assay = "RNA", layer_prefix = "^counts") {
  layers <- grep(layer_prefix, SeuratObject::Layers(seurat_obj[[assay]]), value = TRUE)
  
  pct_list <- lapply(layers, function(ly) {
    M <- SeuratObject::GetAssayData(seurat_obj, assay = assay, layer = ly)
    present <- intersect(target_genes, rownames(M))
    tot <- Matrix::colSums(M)
    tgt <- if (length(present)) Matrix::colSums(M[present, , drop = FALSE]) else rep(0, length(tot))
    setNames(ifelse(tot > 0, 100 * tgt / tot, 0), colnames(M))
  })
  
  pct_vec <- do.call(c, pct_list)
  pct_vec <- pct_vec[colnames(seurat_obj)]
  pct_vec
}


# ---- More helper functions for ggplot ----

make_log_breaks <- function(minv, maxv, include = numeric()) {
  # Add a selected number to the breaks
  # keep only positive values for log scales
  minv <- minv[minv > 0]
  maxv <- maxv[maxv > 0]
  if (!is.finite(minv) || !is.finite(maxv)) return(include)
  
  b <- scales::log_breaks(n = 6)(c(minv, maxv))
  b <- sort(unique(c(b, include)))
  b <- b[b >= minv & b <= maxv]
  b
}

nice_log_labels <- function(x) {
  vapply(x, function(v) {
    if (!is.finite(v)) return(NA_character_)
    if (v == 0) return("0")
    if (abs(v) >= 1) {
      formatC(v, format = "f", digits = 0)          # 100 -> "100"
    } else {
      format(signif(v, 2), scientific = FALSE, trim = TRUE)  # 0.1, 0.01 stay distinct
    }
  }, character(1))
}