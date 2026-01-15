# ---------------------------------------------------------------------
# compute_taxon_contributions
# ---------------------------------------------------------------------
#' @title Compute per-taxon contributions to enzyme load per sample
#' @description
#' Multiply taxon-level enzyme load by taxon abundance per sample to compute
#' a taxon × sample contribution matrix. Optionally normalize abundances to
#' relative proportions (Total Sum Scaling) before multiplication.
#'
#' @param taxon_enzyme_matrix Numeric matrix rows = taxa, cols = enzymes.
#' @param abundance_df Data.frame or matrix rows = taxa, cols = samples.
#' @param enzymes_for_plastic Character vector of enzyme IDs to include or 'All'.
#' @param method_counting_enzyme Character, one of 'abundance' (default) or 'binary'.
#'   - 'abundance': use number of enzyme targets (load) as multiplier of abundance.
#'   - 'binary': treat each taxon-enzyme presence as 0/1 when computing load.
#' @param normalize_taxa_abundance Logical (FALSE) to normalize abundance by simple division by the total count. Default (FALSE).
#'   each sample column sums to 1 (relative contributions). If 'FALSE' use
#'   abundances as provided (absolute contributions).
#'   Skip normalization (set FALSE) when counts are already normalized.
#'
#' @return A list with elements:
#' * 'contributions' : numeric matrix rows = taxa, cols = samples (taxon × sample contribution)
#' * 'enzyme_load'   : numeric vector (named) enzyme load per taxon
#'
#' @details
#' When 'normalize_abundance = TRUE', abundances are converted to proportions
#' via 'A <- sweep(A, 2, colSums(A), "/")'. This corresponds to library-size
#' normalization. Use 'TRUE' when you want comparisons of relative
#' contributions across samples; use 'FALSE' to preserve absolute counts.
#' Skip normalization (set FALSE) when counts are already normalized.
#'
#' @examples
#' # contrib_res <- compute_taxon_contributions(taxon_enzyme_mat, abundance, method_counting_enzyme = "abundance", normalize_taxa_abundance = FALSE)
#'
#' @export
compute_taxon_contributions <- function(taxon_enzyme_matrix,
                                        abundance_df,
                                        enzymes_for_plastic = "All",
                                        method_counting_enzyme = c("abundance", "binary"),
                                        normalize_taxa_abundance = FALSE) {
  
  method <- match.arg(method_counting_enzyme)
  
  # matrix conversion for arithmetic
  A <- as.matrix(abundance_df)
  
  # align taxa: only taxa present in both inputs are used
  common_taxa <- intersect(rownames(taxon_enzyme_matrix), rownames(A))
  if (length(common_taxa) == 0)
    stop("No overlapping taxa between taxon_enzyme_matrix and abundance_df")
  
  # subset and keep structure
  M <- taxon_enzyme_matrix[common_taxa, , drop = FALSE]
  A <- A[common_taxa, , drop = FALSE]
  
  # optional enzyme selection (restrict to a curated set)
  if (!identical(enzymes_for_plastic, "All")) {
    keep <- intersect(colnames(M), enzymes_for_plastic)
    if (length(keep) == 0)
      stop("None of the specified enzymes_for_plastic are present in the enzyme matrix")
    M <- M[, keep, drop = FALSE]
  }
  
  # binary mode: treat >0 as presence
  if (method == "binary") {
    M <- (M > 0) * 1
  }
  
  # enzyme_load per taxon (numeric vector)
  enzyme_load <- rowSums(M)
  
  
  # optionally normalize by sample library size (Total Sum Scaling)
  if (normalize_taxa_abundance) {
    A <- sweep(A, 2, colSums(as.matrix(abundance_df)), "/")
  }
  
  
  # contribution: taxon abundance * enzyme load (per-taxon's contribution to the score)
  contrib <- sweep(A, 1, enzyme_load, "*")
  
  list(contributions = contrib, enzyme_load = enzyme_load)
}
