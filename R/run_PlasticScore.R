# ---------------------------------------------------------------------
# run_pipeline (wrapper)
# ---------------------------------------------------------------------
#' @title Run PlasticScore pipeline
#' @description
#' Wrapper that executes the workflow:
#' parse HMMER outputs, filter hits, build taxon× enzyme matrix, compute plastic
#' scores, compute taxon contributions, run per-taxon tests and prepare objects
#' for plotting.
#'
#' @param hmmsearch_files Character(1). Directory containing '.tbl' / '.domtbl' files.
#' @param abundance_csv Character(1). Path to a CSV file with rows = taxa (rownames), cols = samples.
#' @param metadata_groups_csv Character(1). CSV containing sample metadata (first column = sample IDs).
#' @param metadata_taxa_csv Character(1). CSV mapping 'taxon_id' to taxonomic fields (e.g., Phylum).
#' @param sample_col Optional character(1). Column name in 'metadata_groups_csv' with sample IDs (default = first column).
#' @param condition_col Optional character(1). Column name in 'metadata_groups_csv' with group labels (default = second column).
#' @param min_hmm_coverage Numeric. See 'filter_hmm_hits'.
#' @param max_dom_bias Numeric. See 'filter_hmm_hits'.
#' @param enzymes_for_plastic Character vector or "All". Subset of enzymes to use when scoring.
#' @param method Character. '"abundance"' or '"binary"'.
#' @param normalize_abundance Logical. Passed to 'compute_taxon_contributions' (TSS if TRUE).
#'
#' @return A list (class 'plasticScoreResult') containing:
#' * 'hits_filt', 'taxon_enzyme_mat', 'plastic_score', 'plastic_score_test',
#'   'enzyme_load', 'taxon_stats', 'abundance', 'metadata_groups', 'metadata_taxa',
#'   and plot-ready objects 'df_phylum_enzyme', 'phylum_enzyme_mat', 'enzyme_by_sample',
#'   'contrib_phylum_mat'.
#'
#' @details
#' This wrapper assumes standard input layouts:
#' - 'abundance_csv' should have taxa as rownames and sample IDs as column names.
#' - 'metadata_groups_csv' first two columns are sample ID and group label.
#' - 'metadata_taxa_csv' should map 'taxon_id' to taxonomic ranks such as 'Phylum'.
#'
#' Use 'normalize_abundance = TRUE' to get taxon contributions expressed as
#' proportions (TSS). Use 'FALSE' to preserve absolute abundances.
#'
#' @examples
#' \dontrun{
#' result <- run_pipeline("hmmsearch", "abundance.csv", "metadata_samples.csv", "metadata_taxa.csv")
#' }
#'
#' @export
run_PlasticScore <- function(hmmsearch_files,
                         abundance_csv,
                         metadata_groups_csv,
                         metadata_taxa_csv,
                         sample_col = NULL,
                         condition_col = NULL,
                         min_hmm_coverage = 0.7,
                         max_dom_bias = 0.1,
                         enzymes_for_plastic = "All",
                         method = c("abundance", "binary"),
                         normalize_abundance = TRUE) # relative / compositional (TRUE)
{

  # gather domtbl/tbl files
  hmmtbl_files <- list.files(hmmsearch_files, pattern = "\\.(domtbl|tbl)$", full.names = TRUE)

  # read inputs (abundance rownames must be taxa)
  abundance <- read.csv(abundance_csv, row.names = 1, check.names = FALSE)
  metadata_groups <- read.csv(metadata_groups_csv, stringsAsFactors = FALSE)
  metadata_taxa <- read.csv(metadata_taxa_csv, stringsAsFactors = FALSE)

  if (is.null(sample_col)) sample_col <- colnames(metadata_groups)[1]
  if (is.null(condition_col)) condition_col <- colnames(metadata_groups)[2]

  # parse and filter HMM hits
  hits_raw <- parse_hmmsearch_tbl(hmmtbl_files)
  hits_filt <- filter_hmm_hits(
    hits_raw,
    min_hmm_coverage = min_hmm_coverage,
    max_dom_bias = max_dom_bias
  )

  taxon_enzyme_mat <- build_taxon_enzyme_matrix(hits_filt)

  # compute plastic score and run test
  plastic_score <- compute_plastic_score(
    taxon_enzyme_mat,
    abundance,
    enzymes_for_plastic = enzymes_for_plastic,
    method = method
  )

  plastic_score_test <- compute_plastic_score_test(
    plastic_score, metadata_groups, sample_col, condition_col
  )

  # contributions and per-taxon statistics
  contrib_res <- compute_taxon_contributions(
    taxon_enzyme_mat,
    abundance,
    normalize_abundance = normalize_abundance
  )
  taxon_stats <- test_taxon_contributions(
    contrib_res$contributions, metadata_groups, sample_col, condition_col
  )

  # prepare plot-ready objects: df_phylum_enzyme and phylum_enzyme_mat
  df_phylum_enzyme <- merge(
    hits_filt[, c("taxon_id", "enzyme_id")],
    metadata_taxa,
    by = "taxon_id"
  )
  df_phylum_enzyme <- aggregate(
    taxon_id ~ Phylum + enzyme_id,
    data = df_phylum_enzyme,
    FUN = length
  )
  colnames(df_phylum_enzyme)[3] <- "hit_count"

  phylum_enzyme_mat <- reshape(
    df_phylum_enzyme,
    idvar = "Phylum",
    timevar = "enzyme_id",
    direction = "wide"
  )
  rownames(phylum_enzyme_mat) <- phylum_enzyme_mat$Phylum
  phylum_enzyme_mat <- as.matrix(phylum_enzyme_mat[, -1])
  colnames(phylum_enzyme_mat) <- sub("hit_count.", "", colnames(phylum_enzyme_mat))

  # align taxa and compute enzyme_by_sample via matrix multiplication
  common_taxa <- intersect(
    rownames(taxon_enzyme_mat),
    rownames(abundance)
  )

  if (length(common_taxa) == 0) {
    stop("No overlapping taxa between taxon_enzyme_mat and abundance")
  }

  enzyme_by_sample <- t(
    taxon_enzyme_mat[common_taxa, , drop = FALSE]
  ) %*% as.matrix(
    abundance[common_taxa, , drop = FALSE]
  )

  # aggregate contributions at phylum level for plotting stacked bars
  contrib_phylum_mat <- merge(
    as.data.frame(contrib_res$contributions),
    metadata_taxa,
    by.x = "row.names",
    by.y = "taxon_id"
  )
  contrib_phylum_mat <- aggregate(
    . ~ Phylum,
    data = contrib_phylum_mat[, -1],
    FUN = sum
  )
  rownames(contrib_phylum_mat) <- contrib_phylum_mat$Phylum
  contrib_phylum_mat <- as.matrix(contrib_phylum_mat[, -1])

  # collect result list
  result <- list(
    hits_filt = hits_filt,
    taxon_enzyme_mat = taxon_enzyme_mat,
    plastic_score = plastic_score,
    plastic_score_test = plastic_score_test,
    enzyme_load = contrib_res$enzyme_load,
    taxon_stats = taxon_stats,
    abundance = abundance,
    metadata_groups = metadata_groups,
    metadata_taxa = metadata_taxa,
    df_phylum_enzyme = df_phylum_enzyme,
    phylum_enzyme_mat = phylum_enzyme_mat,
    enzyme_by_sample = enzyme_by_sample,
    contrib_phylum_mat = contrib_phylum_mat
  )

  class(result) <- "plasticScoreResult"
  result
}
