# ---------------------------------------------------------------------
# compute_plastic_score_test
# ---------------------------------------------------------------------
#' @title Wilcoxon test for plastic scores across two groups
#' @description
#' Perform a two-group Wilcoxon (Mann-Whitney) test comparing the plastic score
#' distribution between groups defined in metadata.
#'
#' @param plastic_score Named numeric vector of sample-level scores (names = sample IDs).
#' @param metadata Data.frame containing sample metadata, (columns: Samples and	Groups).
#' @param sample_col Character(1). Column name in 'metadata' that contains sample IDs.
#' @param condition_col Character(1). Column name in 'metadata' that contains group labels.
#'
#' @return An object of class 'htest' produced by 'wilcox.test'.
#' @details
#' The function merges 'plastic_score' with 'metadata' using 'sample_col', then
#' runs 'wilcox.test(score ~ group)'. It is intended for two-group comparisons.
#'
#' @examples
#' # res <- compute_plastic_score_test(score_vector, metadata, "Samples", "Groups")
#'
#' @export
compute_plastic_score_test <- function(plastic_score,
                                       metadata,
                                       sample_col,
                                       condition_col) {

  # assemble data.frame mapping sample -> score
  df <- data.frame(
    sample = names(plastic_score),
    score = as.numeric(plastic_score),
    stringsAsFactors = FALSE
  )

  # merge in metadata to get group labels
  df <- merge(df, metadata, by.x = "sample", by.y = sample_col)

  # create a simple column with group and run Wilcoxon
  df$.group <- df[[condition_col]]

  wilcox.test(score ~ .group, data = df, exact = FALSE)
}
