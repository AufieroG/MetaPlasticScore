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
#' @param normalize_abundance Logical. If 'TRUE' perform column-wise TSS such that
#'   each sample column sums to 1 (relative contributions). If 'FALSE' use
#'   abundances as provided (absolute contributions).
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
#'
#' @examples
#' # contrib_res <- compute_taxon_contributions(taxon_enzyme_mat, abundance, TRUE)
#'
#' @export
compute_taxon_contributions <- function(taxon_enzyme_matrix,
                                        abundance_df,
                                        normalize_abundance = TRUE) {
  
  A <- as.matrix(abundance_df)
  
  # align taxa
  common_taxa <- intersect(rownames(taxon_enzyme_matrix), rownames(A))
  
  # enzyme_load per taxon: sum over all enzymes for that taxon
  enzyme_load <- rowSums(taxon_enzyme_matrix[common_taxa, , drop = FALSE])
  
  # subset abundance to the same taxa
  A <- A[common_taxa, , drop = FALSE]
  
  # optionally normalize by sample library size (Total Sum Scaling)
  if (normalize_abundance) {
    A <- sweep(A, 2, colSums(A), "/")
  }
  
  # contribution: taxon abundance * enzyme load (per-taxon's contribution to the score)
  contrib <- sweep(A, 1, enzyme_load, "*")
  
  list(contributions = contrib, enzyme_load = enzyme_load)
}
