# ---------------------------------------------------------------------
# parse_hmm_filename
# ---------------------------------------------------------------------
#' @title Parse hmmsearch output filename to taxon and enzyme identifiers
#' @description
#' Parse a filename that follows the naming convention produced by the
#' pipeline: '<taxon_id>__<enzyme_id>.tbl' or '<taxon_id>__<enzyme_id>.domtbl'.
#' Returns a list with 'taxon_id' and 'enzyme_id.
#'
#' @param filename Character(1). Path or basename of a hmmsearch output file.
#'
#' @return A named list with elements 'taxon_id' and 'enzyme_id'.
#' @details
#' This helper extracts the taxonomic identifier and enzyme identifier from the
#' filename. It is strict: if the filename does not match the expected pattern an
#' error is raised. 
#'
#' @examples
#' parse_hmm_filename("Bacteroides__PETase.tbl")
#'
#' @export
parse_hmm_filename <- function(filename) {
  # validate input type and length
  if (!is.character(filename) || length(filename) != 1)
    stop("parse_hmm_filename: filename must be a single string")
  
  # take basename and strip recognized extensions
  base <- basename(filename)
  base <- sub("\\.(domtbl|tbl)$", "", base, perl = TRUE)
  
  # split on double underscore to separate taxon and enzyme
  parts <- strsplit(base, "__", fixed = TRUE)[[1]]
  if (length(parts) != 2)
    stop("Filename must follow <taxon_id>__<enzyme_id>.tbl or .domtbl")
  
  # return list for programmatic usage
  list(taxon_id = parts[1], enzyme_id = parts[2])
}
