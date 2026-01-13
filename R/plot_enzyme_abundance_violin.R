#' @title Violin plot of enzyme abundance by condition with Wilcoxon test
#'
#' @description
#' Generate a violin plot (with overlaid boxplot and jittered points) showing
#' raw enzyme abundance across experimental conditions. For each enzyme,
#' a two-group Wilcoxon (Mann–Whitney) test is performed and the FDR-adjusted
#' p-value is annotated on the plot.
#'
#' @details
#' The function expects an enzyme-by-sample abundance matrix (rows = enzymes,
#' columns = samples) and a metadata table linking sample IDs to experimental
#' conditions. The Wilcoxon test is computed independently for each enzyme.
#'
#' Assumptions:
#' \itemize{
#'   \item Exactly two groups are present in \code{condition_col}.
#'   \item Column names of \code{enzyme_by_sample} match sample IDs in \code{meta}.
#' }
#'
#' P-values are adjusted across enzymes using the Benjamini–Hochberg (FDR) method.
#'
#' @param enzyme_by_sample Numeric matrix or data.frame with rows = enzymes and
#'   columns = samples.
#' @param meta Data.frame containing sample metadata.
#' @param sample_col Character(1). Column name in \code{meta} containing sample IDs.
#' @param condition_col Character(1). Column name in \code{meta} defining groups
#'   (two levels expected).
#' @param cfg Optional configuration list (e.g. from \code{default_plot_config()}).
#'
#' @return
#' A ggplot object showing enzyme abundance distributions by condition with
#' Wilcoxon test annotations.
#'
#' @examples
#' \dontrun{
#' cfg <- default_plot_config()
#' p <- plot_enzyme_abundance_violin(
#'   enzyme_by_sample,
#'   metadata_groups,
#'   sample_col = "Samples",
#'   condition_col = "Groups",
#'   cfg = cfg
#' )
#' print(p)
#' }
#'
#' @export
plot_enzyme_abundance_violin <- function(
    enzyme_by_sample,
    meta,
    sample_col,
    condition_col,
    cfg = NULL
) {
  
  # ----------------------------
  # Input validation
  # ----------------------------
  if (!is.matrix(enzyme_by_sample) && !is.data.frame(enzyme_by_sample)) {
    stop("enzyme_by_sample must be a matrix or data.frame (enzymes x samples).")
  }
  if (!is.data.frame(meta)) stop("meta must be a data.frame.")
  if (!sample_col %in% colnames(meta)) stop("sample_col not found in meta.")
  if (!condition_col %in% colnames(meta)) stop("condition_col not found in meta.")
  if (is.null(cfg) || !is.list(cfg)) cfg <- default_plot_config()
  
  # ----------------------------
  # Convert enzyme_by_sample to long format
  # ----------------------------
  df_enzyme <- as.data.frame(enzyme_by_sample)
  df_enzyme$enzyme_id <- rownames(df_enzyme)
  
  df_long <- tidyr::pivot_longer(
    df_enzyme,
    cols = -enzyme_id,
    names_to = "sample",
    values_to = "abundance"
  )
  
  # ----------------------------
  # Merge with metadata
  # ----------------------------
  df_long <- dplyr::left_join(
    df_long,
    meta,
    by = c("sample" = sample_col)
  )
  
  # Drop samples without group assignment
  df_long <- dplyr::filter(df_long, !is.na(.data[[condition_col]]))
  df_long[[condition_col]] <- as.factor(df_long[[condition_col]])
  
  # ----------------------------
  # Wilcoxon test per enzyme
  # ----------------------------
  wilcox_df <- df_long %>%
    dplyr::group_by(enzyme_id) %>%
    dplyr::summarise(
      p_value = stats::wilcox.test(
        abundance ~ .data[[condition_col]],
        exact = FALSE
      )$p.value,
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      p_adj = stats::p.adjust(p_value, method = "fdr"),
      label = paste0("FDR p = ", signif(p_adj, 3))
    )
  
  # Y-position for p-value labels (per enzyme)
  y_pos <- df_long %>%
    dplyr::group_by(enzyme_id) %>%
    dplyr::summarise(
      y = max(abundance, na.rm = TRUE) * 1.05,
      .groups = "drop"
    )
  
  wilcox_df <- dplyr::left_join(wilcox_df, y_pos, by = "enzyme_id")
  
  # ----------------------------
  # Build plot
  # ----------------------------
  p <- ggplot2::ggplot(
    df_long,
    ggplot2::aes(
      x = .data[[condition_col]],
      y = abundance,
      fill = .data[[condition_col]]
    )
  ) +
    ggplot2::geom_violin(trim = TRUE, alpha = 0.5) +
    ggplot2::geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9) +
    ggplot2::geom_jitter(
      width = cfg$jitter_width,
      size = 1.5,
      alpha = cfg$alpha_points
    ) +
    ggplot2::facet_wrap(~ enzyme_id, scales = "free_y") +
    ggplot2::geom_text(
      data = wilcox_df,
      ggplot2::aes(x = 1.5, y = y, label = label),
      inherit.aes = FALSE,
      size = 3
    ) +
    ggplot2::scale_fill_viridis_d(option = cfg$palette_group) +
    ggplot2::theme_minimal(base_size = cfg$base_font_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = cfg$title_size
      ),
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1,
        size = cfg$axis_text_size
      ),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "none"
    ) +
    ggplot2::labs(
      title = "Enzyme abundance by condition",
      x = NULL,
      y = "Abundance"
    )
  
  p
}
