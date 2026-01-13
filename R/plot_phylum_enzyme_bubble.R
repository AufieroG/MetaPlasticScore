#' @title Bubble plot of enzyme hits by phylum
#'
#' @description
#' Create a bubble plot where each point represents the number of HMM hits
#' for a given (enzyme, phylum) pair. Bubble area encodes hit_count and color
#' encodes phylum.
#'
#' @param df_phylum_enzyme Data.frame with columns: Phylum, enzyme_id, hit_count.
#' @param cfg Configuration list produced by \code{default_plot_config()}.
#'
#' @return ggplot object (bubble plot)
#'
#' @examples
#' \dontrun{
#' cfg <- default_plot_config()
#' plot_phylum_enzyme_bubble(df_phylum_enzyme, cfg)
#' }
#'
#' @export
plot_phylum_enzyme_bubble <- function(df_phylum_enzyme, cfg) {
  # validate input type - fail early with informative message
  if (!is.data.frame(df_phylum_enzyme)) {
    stop("df_phylum_enzyme must be a data.frame with columns Phylum, enzyme_id, hit_count")
  }
  # ensure cfg exists and has defaults if missing
  if (missing(cfg) || !is.list(cfg)) cfg <- default_plot_config()
  
  # Build ggplot mapping:
  # x = enzyme identifier, y = Phylum (categorical),
  # size = hit_count (bubble area), color = Phylum
  p <- ggplot2::ggplot(
    df_phylum_enzyme,
    ggplot2::aes(x = enzyme_id, y = Phylum, size = hit_count, color = Phylum)
  ) +
    # geom_point with alpha to reduce overplotting; size aesthetic maps to area by default
    ggplot2::geom_point(alpha = cfg$alpha_points) +
    # color scale for discrete Phylum; use Viridis discrete for accessibility
    ggplot2::scale_color_viridis_d(option = cfg$palette_group) +
    # scale sizes so very large counts don't produce giant bubbles
    ggplot2::scale_size_area(max_size = 12) +
    # base minimal theme for clean appearance
    ggplot2::theme_minimal(base_size = cfg$base_font_size) +
    # theme tweaks: center title, rotate x labels (likely many enzymes), size axes labels
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = cfg$title_size),
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = cfg$axis_text_size),
      axis.text.y = ggplot2::element_text(size = cfg$axis_text_size)
    ) +
    # labels
    ggplot2::labs(title = "Bubble plot: enzyme hits per Phylum", x = "Enzyme", y = NULL)
  
  # Return ggplot object for further customization by caller
  p
}