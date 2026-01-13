#' @title PCoA of enzyme functional profiles (PCoA via cmdscale)
#'
#' @description
#' Calculate and plot a 2D PCoA (classical MDS) of samples using enzyme-by-sample
#' profiles. The function supports optional transformations (none/log/hellinger)
#' and two distance options (Bray or Euclidean). Sample grouping (color/fill)
#' is taken from the provided metadata (first column = sample id, second column = group).
#'
#' @param enzyme_by_sample Numeric matrix or data.frame with rows = enzymes, cols = samples.
#' @param meta Data.frame with sample metadata. First column must contain sample IDs
#'             matching the column names of \code{enzyme_by_sample}. Second column is used as group.
#' @param cfg Configuration list produced by \code{default_plot_config()}.
#'
#' @return ggplot object (PCoA scatter plot)
#'
#' @examples
#' \dontrun{
#' cfg <- default_plot_config()
#' cfg$pcoa_transform <- "hellinger"
#' cfg$pcoa_distance <- "bray"
#' plot_pcoa_enzyme_profiles(enzyme_by_sample, metadata_groups, cfg)
#' }
#'
#' @export
plot_pcoa_enzyme_profiles <- function(enzyme_by_sample, meta, cfg) {
  # require vegan package for decostand and vegdist (if used)
  if (!requireNamespace("vegan", quietly = TRUE)) {
    stop("Package 'vegan' is required for transformations and Bray distance. Install it before running.")
  }
  # ensure configuration is available
  if (missing(cfg) || !is.list(cfg)) cfg <- default_plot_config()
  
  # convert input to matrix and check orientation
  mat <- as.matrix(enzyme_by_sample)
  # Expect enzymes x samples; we need samples × enzymes for distance calculations
  if (ncol(mat) < 2) {
    stop("Need at least two samples (columns) in enzyme_by_sample to compute PCoA.")
  }
  
  # read transform & distance parameters from cfg (with safe defaults)
  transform <- if (!is.null(cfg$pcoa_transform)) cfg$pcoa_transform else "none"
  distance  <- if (!is.null(cfg$pcoa_distance)) cfg$pcoa_distance else "bray"
  
  # transpose matrix to samples x enzymes (rows = samples)
  samp_mat <- t(mat)
  
  # ----------------------------
  # Apply requested transformation
  # ----------------------------
  # 'none' : raw counts (use only if comparable across samples / already normalized)
  # 'log'  : log1p reduces influence of very abundant enzymes
  # 'hellinger' : sqrt of relative abundances; recommended before Euclidean for compositional data
  if (transform == "log") {
    # log1p is robust to zeros and widely used for count-like data prior to distance calc
    samp_mat <- log1p(samp_mat)
  } else if (transform == "hellinger") {
    # convert to relative abundances per sample and then take sqrt; vegan::decostand handles this
    samp_mat <- vegan::decostand(samp_mat, method = "hellinger")
  } else if (transform != "none") {
    stop("cfg$pcoa_transform must be one of: 'none', 'log', 'hellinger'")
  }
  
  # ----------------------------
  # Compute distance matrix
  # ----------------------------
  # Bray is default and widely used for ecological/functional profiles; Euclidean supported for transformed data
  if (distance == "bray") {
    dist_samp <- vegan::vegdist(samp_mat, method = "bray")
  } else if (distance == "euclidean") {
    dist_samp <- stats::dist(samp_mat) # base R Euclidean distance
  } else {
    stop("cfg$pcoa_distance must be one of: 'bray', 'euclidean'")
  }
  
  # ----------------------------
  # Classical MDS / PCoA
  # ----------------------------
  # stats::cmdscale returns coordinates and eigenvalues
  pcoa <- stats::cmdscale(dist_samp, k = 2, eig = TRUE)
  
  # Extract eigenvalues to compute percent variance explained (positive eigenvalues only)
  eig <- pcoa$eig
  # defensive: if eig is NULL or contains NA, set to zeros to avoid errors
  if (is.null(eig) || all(is.na(eig))) {
    warning("PCoA returned no eigenvalues; variance percentages set to NA")
    var1 <- NA
    var2 <- NA
  } else {
    pos_sum <- sum(eig[eig > 0])
    # protect against division by zero
    if (pos_sum <= 0) {
      var1 <- NA
      var2 <- NA
    } else {
      var1 <- round(100 * eig[1] / pos_sum, 1)
      var2 <- round(100 * eig[2] / pos_sum, 1)
    }
  }
  
  # ----------------------------
  # Prepare plotting data.frame by merging metadata
  # ----------------------------
  ord_df <- data.frame(
    Sample = rownames(pcoa$points),
    PCoA1 = pcoa$points[, 1],
    PCoA2 = pcoa$points[, 2],
    stringsAsFactors = FALSE
  )
  
  # If meta provided, merge to attach group information (assume meta first column = sample id)
  if (!is.null(meta) && is.data.frame(meta) && ncol(meta) >= 2) {
    plot_df <- merge(ord_df, meta, by.x = "Sample", by.y = colnames(meta)[1], all.x = TRUE)
    group_col <- colnames(meta)[2]
    # ensure group is a factor for controlled plotting
    plot_df[[group_col]] <- as.factor(plot_df[[group_col]])
  } else {
    # fallback: create empty group column so aesthetics function later does not fail
    plot_df <- ord_df
    group_col <- NULL
  }
  
  # ----------------------------
  # Plot
  # ----------------------------
  # Use shape 21 (point with fill + stroke) so we can use fill=group and black outline
  if (!is.null(group_col)) {
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = PCoA1, y = PCoA2, fill = .data[[group_col]])) +
      ggplot2::geom_point(shape = 21, color = "black", stroke = 0.3, size = cfg$pcoa_point_size, alpha = cfg$alpha_points) +
      ggplot2::scale_fill_viridis_d(option = cfg$palette_group, name = group_col)
  } else {
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = PCoA1, y = PCoA2)) +
      ggplot2::geom_point(size = cfg$pcoa_point_size, alpha = cfg$alpha_points)
  }
  
  # theme and labels (center title, show % variance in axis labels when available)
  p <- p +
    ggplot2::theme_minimal(base_size = cfg$base_font_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = cfg$title_size),
      axis.text = ggplot2::element_text(size = cfg$axis_text_size)
    ) +
    ggplot2::labs(
      title = "PCoA of samples by enzyme functional profile",
      subtitle = paste0("Distance: ", distance, "; Transform: ", transform),
      x = ifelse(is.na(var1), "PCoA1", paste0("PCoA1 (", var1, "%)")),
      y = ifelse(is.na(var2), "PCoA2", paste0("PCoA2 (", var2, "%)"))
    )
  
  # add labels if requested and ggrepel available
  if (isTRUE(cfg$pcoa_label) && requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = Sample), size = cfg$axis_text_size / 3)
  } else if (isTRUE(cfg$pcoa_label)) {
    # fallback to geom_text if ggrepel missing (might overlap)
    p <- p + ggplot2::geom_text(ggplot2::aes(label = Sample), vjust = -1.2, size = cfg$axis_text_size / 3)
  }
  
  # Return ggplot object
  p
}