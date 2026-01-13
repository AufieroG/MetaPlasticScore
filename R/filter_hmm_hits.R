# ---------------------------------------------------------------------
# filter_hmm_hits
# ---------------------------------------------------------------------
#' @title Filter HMM hits by coverage and domain bias
#' @description
#' Apply quality filters to the parsed hmmsearch hit table. The default
#' filters keep hits with HMM coverage >= 0.9 and domain bias <= 0.1.
#'
#' @param hits Data.frame returned by 'parse_hmmsearch_tbl'.
#' @param min_hmm_coverage Numeric scalar. Minimum fraction of the HMM that must be aligned.
#' @param max_dom_bias Numeric scalar. Maximum allowed domain bias (NA allowed).
#'
#' @return Subset of 'hits' meeting the filter criteria.
#' @details
#' - 'min_hmm_coverage' helps exclude partial / truncated HMM matches
#' - 'max_dom_bias' is intended to remove hits where the HMM scoring is biased
#'   by low-complexity regions (higher bias values are more suspicious)
#'
#'
#' @examples
#' # hits <- parse_hmmsearch_tbl(files)
#' # hits_filt <- filter_hmm_hits(hits, min_hmm_coverage = 0.8, max_dom_bias = 0.05)
#'
#' @export
filter_hmm_hits <- function(hits, min_hmm_coverage = 0.9, max_dom_bias = 0.1) {
  subset(
    hits,
    !is.na(hmm_coverage) &
      hmm_coverage >= min_hmm_coverage &
      (is.na(dom_bias) | dom_bias <= max_dom_bias)
  )
}