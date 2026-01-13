#' @title Default plotting configuration
#'
#' @description
#' Centralized configuration controlling appearance and behavior of all plotting
#' functions in the Plastic Score pipeline.
#'
#' This function returns a named list of parameters that can be modified by the
#' user and passed unchanged to all plotting functions, ensuring visual and
#' methodological consistency across figures.
#'
#' Options are grouped by plot type to make their purpose explicit.
#'
#' @return
#' A named list of plotting parameters.
#'
#' @export
default_plot_config <- function() {
  
  list(
    # ============================================================
    # OUTPUT / SAVING OPTIONS (GLOBAL)
    
    output_dir = "Plots",      # Directory where plots are saved if save = TRUE
    save = TRUE,               # If TRUE, plots are written to disk automatically
    dpi = 300,                 # Resolution for raster images (publication-ready)
    width = 8,                 # Default plot width in inches
    height = 5,                # Default plot height in inches
    
    # ============================================================
    # GENERAL SIZES
    
    base_font_size = 10,       # Base font size used by theme_minimal()
    axis_text_size = 9,       # Font size for axis tick labels
    title_size = 14,           # Font size for plot titles
    
    # ============================================================
    # COLOR PALETTES
    
    palette_enzyme = "cividis", # Used for enzyme-level stacked bars
    palette_group  = "cividis", # Used for group-based plots (PCoA, boxplots)
    palette_phylum = "Set2",    # Used for phylum-level stacked contributions
    
    # ============================================================
    # GENERAL AESTHETICS
    
    alpha_points = 0.8,        # Transparency for scatter points
    jitter_width = 0.15,       # Horizontal jitter for overplotted points
    
    # ============================================================
    # HEATMAP OPTIONS (plot_enzyme_heatmap)
    
    heatmap_row_level = "Taxon",   # "Taxon" or "Phylum"
    heatmap_log = TRUE,            # Apply log10(x + 1) for visualization only
    heatmap_cluster_rows = TRUE,   # Hierarchical clustering of rows
    heatmap_cluster_cols = TRUE,   # Hierarchical clustering of columns
    heatmap_top_n_taxa = NULL,     # Show only top-N rows (by total abundance)
    heatmap_title = "Enzyme abundance heatmap",
    heatmap_row_label_size = 8,    # Font size of row labels
    heatmap_col_label_size = 8,    # Font size of column labels
    
    # ============================================================
    # PCoA OPTIONS (plot_pcoa_enzyme_profiles)
    
    pcoa_label = TRUE,             # Draw sample labels using ggrepel
    pcoa_point_size = 3,           # Size of points in ordination plot
    pcoa_transform = "none",       # "none" | "log"
    pcoa_distance  = "bray",       # "bray" | "euclidean"
    
    # ============================================================
    # STACKED CONTRIBUTION PLOTS
    
    normalize_contrib = TRUE      # Show contributions as percentages per sample
    
  )

}
