# ---------------------------------------------------------------------
# build_taxon_enzyme_matrix
# ---------------------------------------------------------------------
#' @title Build taxon × enzyme matrix from HMM hits
#' @description
#' Aggregate parsed HMM hits to produce a matrix whose rows are taxon identifiers
#' and columns are enzyme identifiers. Cell values are counts of unique target
#' sequences observed for that taxon-enzyme pair.
#'
#' @param hits Data.frame produced by 'parse_hmmsearch_tbl' (preferably filtered).
#'
#' @return Numeric matrix with rownames = taxon_id and colnames = enzyme_id.
#' @details
#' The function:
#' - deduplicates 'target_name' per taxon-enzyme pair,
#' - counts unique targets per taxon × enzyme,
#' - constructs a dense matrix filled with zeros for missing pairs.
#'
#' Use this matrix to compute per-taxon enzyme load (row sums) or to map enzymes
#' to taxa for downstream calculations.
#'
#' @examples
#' # agg_mat <- build_taxon_enzyme_matrix(hits_filt)
#'
#' @export
build_taxon_enzyme_matrix <- function(hits) {
  
  # keep only relevant columns and remove duplicates
  df <- unique(hits[, c("taxon_id", "enzyme_id", "target_name")])
  
  # aggregate by taxon_id + enzyme_id counting unique target_name
  agg <- aggregate(
    target_name ~ taxon_id + enzyme_id,
    data = df,
    FUN = function(x) length(unique(x))
  )
  
  # unique taxa and enzymes
  taxa <- unique(agg$taxon_id)
  enzymes <- unique(agg$enzyme_id)
  
  # initialize 0 matrix
  mat <- matrix(
    0,
    nrow = length(taxa),
    ncol = length(enzymes),
    dimnames = list(taxa, enzymes)
  )
  
  # populate matrix with aggregated counts
  for (i in seq_len(nrow(agg))) {
    mat[agg$taxon_id[i], agg$enzyme_id[i]] <- agg$target_name[i]
  }
  
  mat
}
