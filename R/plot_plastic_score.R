#' @title Plot plastic score by group (violin + box + jitter)
#'
#' @description
#' Create a combined violin + boxplot + jitter plot showing plastic score per group.
#' The function expects a named numeric vector of scores and metadata linking samples to groups.
#'
#' @param plastic_score Named numeric vector (names = sample IDs).
#' @param meta Data.frame with sample metadata (must contain sample and grouping columns).
#' @param sample_col Character name of column in \code{meta} that matches names(plastic_score).
#' @param condition_col Character name of grouping column in \code{meta}.
#' @param cfg Configuration list produced by \code{default_plot_config()}.
#' @param p_value Optional numeric p-value to annotate the plot.
#'
#' @return ggplot object (violin + boxplot)
#'
#' @examples
#' \dontrun{
#' plot_plastic_score(plastic_score, metadata_groups, "Samples", "Groups", cfg)
#' }
#'
#' @export
plot_plastic_score <- function(plastic_score, meta, sample_col, condition_col, cfg, p_value = NULL) {
  # Input validation: plastic_score must be named numeric vector
  if (!is.numeric(plastic_score) || is.null(names(plastic_score))) {
    stop("plastic_score must be a named numeric vector (names = sample IDs).")
  }
  if (!is.data.frame(meta)) stop("meta must be a data.frame with sample metadata.")
  if (!sample_col %in% colnames(meta)) stop("sample_col not found in meta.")
  if (!condition_col %in% colnames(meta)) stop("condition_col not found in meta.")
  if (missing(cfg) || !is.list(cfg)) cfg <- default_plot_config()
  
  # Merge scores into metadata to prepare plotting table
  df_score <- data.frame(sample = names(plastic_score), score = as.numeric(plastic_score), stringsAsFactors = FALSE)
  plot_df <- merge(df_score, meta, by.x = "sample", by.y = sample_col, all.x = TRUE)
  
  # warn if some samples lack metadata
  if (any(is.na(plot_df[[condition_col]]))) {
    warning("Some samples in plastic_score do not have matching metadata rows; they will appear with NA grouping.")
  }
  
  # build ggplot with aes_string to allow column name in variable
  p <- ggplot2::ggplot(plot_df, ggplot2::aes_string(x = condition_col, y = "score", fill = condition_col)) +
    # violin for distribution (alpha lower so jitter/box are visible)
    ggplot2::geom_violin(alpha = 0.3, trim = TRUE) +
    # narrow boxplot to show median/IQR on top of violin
    ggplot2::geom_boxplot(width = 0.15, outlier.shape = NA) +
    # add jittered points to show individual samples
    ggplot2::geom_jitter(width = cfg$jitter_width, alpha = cfg$alpha_points) +
    # fill palette for groups
    ggplot2::scale_fill_viridis_d(option = cfg$palette_group) +
    ggplot2::theme_minimal(base_size = cfg$base_font_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = cfg$title_size),
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = cfg$axis_text_size),
      axis.text.y = ggplot2::element_text(size = cfg$axis_text_size)
    ) +
    ggplot2::labs(title = "Plastic score by group", y = "Plastic degradation potential", x = NULL)
  
  # annotate p-value near top of plot if provided
  if (!is.null(p_value) && is.numeric(p_value)) {
    # place annotation at ~5% above maximum score (display only)
    y_pos <- max(plot_df$score, na.rm = TRUE) * 1.05
    p <- p + ggplot2::annotate("text", x = 1.5, y = y_pos, label = paste0("p = ", signif(p_value, 3)))
  }
  
  p
}