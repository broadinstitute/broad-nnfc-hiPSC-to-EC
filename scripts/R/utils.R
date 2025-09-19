# ---- Color Palettes ----

#' Qualitative color palette (4 colors)
#'
#' Useful for categorical data visualization.
palette_4_qual <- c("#1B9E77", "#D95F02", "#7570B3", "#E6AB02")

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
read_10x_rna_grna <- function(h5, sample_name) {
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
  if (!is.na(grna_name)) {
    obj[["gRNA"]] <- Seurat::CreateAssayObject(counts = x[[grna_name]])
  } else {
    message("No CRISPR/gRNA modality found in: ", h5)
  }
  
  obj
}