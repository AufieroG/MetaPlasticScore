#' @title Wrapper: generate and optionally save all plots
#'
#' @description
#' Convenience wrapper that calls all plotting functions (where inputs are present),
#' collects produced plot objects and saves ggplot objects to disk when requested.
#' ComplexHeatmap objects are not saved by ggsave and the function prints instructions
#' for saving heatmaps to file.
#'
#' @param df_phylum_enzyme Data.frame for phylum vs enzyme stacked bar / bubble plots.
#' @param taxon_enzyme_mat Matrix taxon x enzyme (rows taxa, cols enzymes).
#' @param phylum_enzyme_mat Matrix phylum x enzyme.
#' @param enzyme_by_sample Matrix enzymes x samples (rows enzymes, cols samples).
#' @param plastic_score Named numeric vector (names = sample IDs).
#' @param contrib_phylum_mat Numeric matrix rows = Phylum, cols = Sample.
#' @param meta Data.frame with sample metadata (first column sample IDs, second column group).
#' @param meta_taxa Data.frame mapping taxon_id -> Phylum.
#' @param sample_col Name of sample ID column in meta (used by plastic score test).
#' @param condition_col Name of grouping column in meta (used by plastic score test).
#' @param cfg Configuration list produced by \code{default_plot_config()}.
#'
#' @return Named list of plot objects (invisible). Heatmap element (if present) is a ComplexHeatmap object.
#'
#' @export
plot_MetaPlasticScore <- function(
    df_phylum_enzyme = NULL,
    taxon_enzyme_mat = NULL,
    phylum_enzyme_mat = NULL,
    enzyme_by_sample = NULL,
    plastic_score = NULL,
    contrib_phylum_mat = NULL,
    meta = NULL,
    meta_taxa = NULL,
    sample_col = NULL,
    condition_col = NULL,
    cfg = NULL
) {
  # ensure cfg present
  if (is.null(cfg) || !is.list(cfg)) cfg <- default_plot_config()
  
  # create output directory if save requested
  if (isTRUE(cfg$save) && !dir.exists(cfg$output_dir)) {
    dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  plots <- list()
  
  # 1) stacked bar - require df_phylum_enzyme
  plots$stacked_bar <- tryCatch({
    if (is.null(df_phylum_enzyme)) stop("df_phylum_enzyme required for stacked bar")
    plot_phylum_enzyme_bar(df_phylum_enzyme, cfg)
  }, error = function(e) {
    message("stacked_bar skipped: ", e$message); NULL
  })
  
  # 2) heatmap - need taxon_enzyme_mat or phylum_enzyme_mat
  plots$heatmap <- tryCatch({
    if (!is.null(taxon_enzyme_mat) || !is.null(phylum_enzyme_mat)) {
      plot_enzyme_heatmap(taxon_enzyme_mat, phylum_enzyme_mat, meta_taxa, cfg, draw = FALSE)
    } else {
      message("heatmap skipped: neither taxon_enzyme_mat nor phylum_enzyme_mat provided."); NULL
    }
  }, error = function(e) {
    message("heatmap skipped: ", e$message); NULL
  })
  
  # 3) bubble
  plots$bubble <- tryCatch({
    if (is.null(df_phylum_enzyme)) stop("df_phylum_enzyme required for bubble plot")
    plot_phylum_enzyme_bubble(df_phylum_enzyme, cfg)
  }, error = function(e) {
    message("bubble skipped: ", e$message); NULL
  })
  
  # 4) PCoA - requires enzyme_by_sample (enzymes x samples)
  plots$pcoa <- tryCatch({
    if (is.null(enzyme_by_sample)) stop("enzyme_by_sample required for PCoA")
    mat <- as.matrix(enzyme_by_sample)
    if (ncol(mat) < 2) stop("PCoA requires at least 2 samples (columns) in enzyme_by_sample.")
    plot_pcoa_enzyme_profiles(mat, meta, cfg)
  }, error = function(e) {
    message("pcoa skipped: ", e$message); NULL
  })
  
  # 5) plastic score
  plots$plastic_score <- tryCatch({
    if (is.null(plastic_score) || is.null(meta) || is.null(sample_col) || is.null(condition_col)) {
      stop("Provide plastic_score, meta, sample_col, and condition_col for plastic score plot")
    }
    plot_plastic_score(plastic_score, meta, sample_col, condition_col, cfg)
  }, error = function(e) {
    message("plastic_score skipped: ", e$message); NULL
  })
  
  # 6) contributions per phylum
  plots$contributions <- tryCatch({
    if (is.null(contrib_phylum_mat)) stop("contrib_phylum_mat required for contributions plot")
    plot_contributions_phylum(contrib_phylum_mat, cfg)
  }, error = function(e) {
    message("contributions skipped: ", e$message); NULL
  })
  
  # # 7) plot_enzyme_abundance_violin
  plots$abundance_violin <- tryCatch({
    if (is.null(enzyme_by_sample)) stop("enzyme_by_sample required for plot_enzyme_abundance_violin plot")
    plot_enzyme_abundance_violin(enzyme_by_sample, meta, sample_col, condition_col, cfg )
  }, error = function(e) {
    message("plot_enzyme_abundance_violin skipped: ", e$message); NULL
  })
  
  # Save ggplot objects (ComplexHeatmap objects handled separately)
  if (isTRUE(cfg$save)) {
    for (nm in names(plots)) {
      obj <- plots[[nm]]
      if (is.null(obj)) next
      # ggplot objects: save as PNG
      if (inherits(obj, "ggplot")) {
        fname <- file.path(cfg$output_dir, paste0(nm, ".png"))
        tryCatch({
          ggplot2::ggsave(filename = fname, plot = obj, width = cfg$width, height = cfg$height, dpi = cfg$dpi)
        }, error = function(e) {
          message("Failed to save plot '", nm, "': ", e$message)
        })
      } else if (inherits(obj, "Heatmap") || inherits(obj, "HeatmapList")) {
        # for ComplexHeatmap, provide explicit instructions (cannot use ggsave)
        message("Plot '", nm, "' is a ComplexHeatmap object. To save it, run:\n",
                "  pdf('", file.path(cfg$output_dir, paste0(nm, '.pdf')), "', width=", cfg$width, ", height=", cfg$height, ")\n",
                "  ComplexHeatmap::draw(plots$heatmap)\n",
                "  dev.off()")
      } else {
        message("Skipping saving for plot '", nm, "' of class ", paste(class(obj), collapse = "/"))
      }
    }
  }
  
  invisible(plots)
}