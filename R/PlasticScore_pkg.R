#' plasticScorePipeline.R
#'
#' Pipeline for computing a "plastic degradation potential"
#' score from HMMER (hmmsearch) outputs and taxon abundance tables.
#'
#' The file contains parsers for HMMER output, matrix builders, scoring logic,
#' statistical tests and plotting helpers (ggplot2 + ComplexHeatmap).
#'
#'
#' @import stats
#' @import utils
#' @import ggplot2
#' @import vegan
#' @import tidyr
#' @import circlize
#' @import grid
#' @import ggrepel
#' @import viridis
#' @import RColorBrewer
#' @import magrittr
#' @import dplyr
#'
#' @author
#' Gaetano Aufiero
#'
#' Maintainer: Gaetano Aufiero
#' @name PlasticScore
NULL
