library("MetaPlasticScore")

test_that("PlasticScore runs on example data in inst/extdata", {
# locate example data shipped with the package
extdata_path <- system.file("extdata", package = "MetaPlasticScore")
expect_true(dir.exists(extdata_path))

hmm_dir <- file.path(extdata_path, "hmmsearch")
abundance_csv <- file.path(extdata_path, "Abundance_simulated.csv")
metadata_groups_csv <- file.path(extdata_path, "Metadata_groups.csv")
metadata_taxa_csv <- file.path(extdata_path, "Metadata_taxa.csv")


# run the pipeline on example data
result <- run_MetaPlasticScore(
  hmmsearch_files = hmm_dir,
  abundance_csv = abundance_csv,
  metadata_groups_csv = metadata_groups_csv,
  metadata_taxa_csv = metadata_taxa_csv,
  min_hmm_coverage = 0.7,
  max_dom_bias = 0.1,
  enzymes_for_plastic = "All",
  method = "abundance",
  normalize_abundance = FALSE
)

# check output class
expect_s3_class(result, "MetaPlasticScoreResult")

})
