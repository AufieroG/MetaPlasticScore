# ---------------------------------------------------------------------
# compute_plastic_score
# ---------------------------------------------------------------------
#' @title Compute sample-level plastic degradation potential score
#' @description
#' Compute a per-sample score representing the potential for plastic degradation
#' given a taxon × enzyme matrix and a taxon abundance table.
#'
#' @param taxon_enzyme_matrix Numeric matrix rows = taxa, cols = enzymes.
#'   Row names must match taxon IDs used in 'abundance_df'.
#' @param abundance_df Data.frame or matrix with rows = taxa, cols = samples.
#'   Row names must be taxon IDs.
#' @param enzymes_for_plastic Character vector of enzyme IDs to include or 'All'.
#' @param method Character, one of 'abundance' (default) or 'binary'.
#'   - 'abundance': use number of enzyme targets (load) as multiplier of abundance.
#'   - 'binary': treat each taxon-enzyme presence as 0/1 when computing load.
#'
#' @return Named numeric vector of scores (names = sample IDs).
#' @details
#' Steps:
#' 1. Intersect taxa between 'taxon_enzyme_matrix' and 'abundance_df'.
#' 2. Optionally subselect enzymes via 'enzymes_for_plastic'.
#' 3. If 'method = binary', convert enzyme counts to presence/absence per taxon.
#' 4. Compute enzyme load per taxon as row sums of the enzyme matrix.
#' 5. Score each sample as sum over taxa of (abundance × enzyme_load).
#'
#' Interpretation notes:
#' - If 'abundance_df' is library-size normalized (e.g., columns sum = 1),
#'   the score is relative. If abundances are absolute read counts, the score is
#'   absolute and depends on library size.
#' - Consider normalizing abundance consistently across samples if you intend to
#'   compare scores between samples.
#'
#' @examples
#' # score <- compute_plastic_score(taxon_enzyme_mat, abundance, method = "abundance")
#'
#' @export
compute_plastic_score <- function(taxon_enzyme_matrix,
                                  abundance_df,
                                  enzymes_for_plastic = "All",
                                  method = c("abundance", "binary")) {
  method <- match.arg(method)
  
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
  
  # multiply abundances by enzyme load and sum by column (sample)
  score <- colSums(A * enzyme_load)
  
  names(score) <- colnames(A)
  score
}