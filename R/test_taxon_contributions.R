# ---------------------------------------------------------------------
# test_taxon_contributions
# ---------------------------------------------------------------------
#' @title Statistical tests on taxon contributions between two groups
#' @description
#' For each taxon, run a two-sample Wilcoxon test comparing its contribution
#' across two groups of samples defined in the provided metadata.
#'
#' @param contribution_matrix Numeric matrix rows = taxa, cols = samples (same ordering as metadata).
#' @param metadata Data.frame containing sample metadata.
#' @param sample_col Character(1). Column name of sample IDs in metadata.
#' @param condition_col Character(1). Column name with group labels (two unique values expected).
#' @param pseudo Numeric. Pseudocount added to means to compute log2 fold changes safely.
#'
#' @return Data.frame with columns: 'taxon', 'mean_group1', 'mean_group2', 'log2FC', 'p_value', 'p_adj'.
#' @details
#' - Samples are aligned by matching column names of 'contribution_matrix' to 'metadata'.
#' - 'lev <- unique(groups)' is used to pick group order for reporting 'mean_group1' and 'mean_group2'.
#' - Taxa with all-zero values in both groups are skipped (returned as NULL).
#' - p-values are adjusted for multiple testing using FDR (Benjamini-Hochberg).
#'
#' @examples
#' # res_taxa <- test_taxon_contributions(contrib_matrix, metadata, "Samples", "Groups")
#'
#' @export
test_taxon_contributions <- function(contribution_matrix,
                                     metadata,
                                     sample_col,
                                     condition_col,
                                     pseudo = 1e-9) {

  samples <- colnames(contribution_matrix)
  taxa <- rownames(contribution_matrix)

  # align metadata row order to sample columns
  meta <- metadata[match(samples, metadata[[sample_col]]), ]
  groups <- meta[[condition_col]]
  lev <- unique(groups)

  # compute per-taxon comparisons
  res <- lapply(taxa, function(tax) {
    x <- contribution_matrix[tax, ]
    g1 <- x[groups == lev[1]]
    g2 <- x[groups == lev[2]]

    if (all(g1 == 0) & all(g2 == 0)) return(NULL)

    wt <- wilcox.test(g1, g2, exact = FALSE)

    data.frame(
      taxon = tax,
      mean_group1 = mean(g1),
      mean_group2 = mean(g2),
      log2FC = log2((mean(g1) + pseudo) / (mean(g2) + pseudo)),
      p_value = wt$p.value,
      stringsAsFactors = FALSE
    )
  })

  res <- do.call(rbind, res)
  res$p_adj <- p.adjust(res$p_value, method = "fdr")
  res
}
