#' @title Stacked bar of phylum contributions per sample
#'
#' @description
#' Plot contributions aggregated at phylum level for each sample. Optionally
#' normalize per sample (column) so contributions are shown as percentages.
#'
#' @param contrib_phylum_mat Numeric matrix with rows = Phylum, cols = Sample.
#' @param cfg Configuration list produced by \code{default_plot_config()}.
#'
#' @return ggplot object (stacked bar)
#'
#' @examples
#' \dontrun{
#' plot_contributions_phylum(contrib_phylum_mat, cfg)
#' }
#'
#' @export
plot_contributions_phylum <- function(contrib_phylum_mat, cfg) {
  # input checks
  if (!is.matrix(contrib_phylum_mat) && !is.data.frame(contrib_phylum_mat)) {
    stop("contrib_phylum_mat must be a matrix or data.frame with rows=Phylum, cols=Sample.")
  }
  if (missing(cfg) || !is.list(cfg)) cfg <- default_plot_config()
  
  # convert to matrix for numeric operations
  mat <- as.matrix(contrib_phylum_mat)
  
  # Normalize per sample to percent if requested: divide each column by its sum
  if (isTRUE(cfg$normalize_contrib)) {
    col_sums <- colSums(mat, na.rm = TRUE)
    # avoid division by zero by replacing zero sums with 1 (resulting percentages = 0)
    col_sums[col_sums == 0] <- 1
    mat <- sweep(mat, 2, col_sums, "/") * 100
  }
  
  # Convert to long format for ggplot: we need Sample, Phylum, value columns
  df <- as.data.frame(t(mat), stringsAsFactors = FALSE) # rows now samples
  df$Sample <- rownames(df)
  # pivot_longer (tidyr) to create long tidy dataframe
  df_long <- tidyr::pivot_longer(df, -Sample, names_to = "Phylum", values_to = "value")
  
  # build stacked bar chart
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = Sample, y = value, fill = Phylum)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_brewer(palette = cfg$palette_phylum) +
    ggplot2::theme_minimal(base_size = cfg$base_font_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = cfg$title_size),
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = cfg$axis_text_size),
      axis.text.y = ggplot2::element_text(size = cfg$axis_text_size)
    ) +
    ggplot2::labs(
      title = "Phylum contributions per sample",
      y = ifelse(cfg$normalize_contrib, "Percent contribution", "Contribution"),
      x = NULL
    )
  
  p
}