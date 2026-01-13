# ---------------------------------------------------------------------
# parse_hmmsearch_tbl
# ---------------------------------------------------------------------
#' @title Parse hmmsearch domtbl/tbl outputs into a data.frame of hits
#' @description
#' Read a set of 'domtbl'/'tbl' files produced by HMMER (hmmsearch) and parse
#' the main tabular fields into a single data.frame. Files are expected to be
#' named following '<taxon_id>__<enzyme_id>.tbl' (or '.domtbl') so taxon and
#' enzyme identifiers can be inferred automatically.
#'
#' @param files Character vector. Paths to '.tbl' or '.domtbl' files.
#'
#' @return A data.frame containing parsed columns for each hit:
#' * 'target_name', 'tlen', 'query_name', 'qlen', 'full_evalue', 'full_score',
#'   'dom_num', 'dom_of', 'i_evalue', 'dom_score', 'dom_bias', 'hmm_from',
#'   'hmm_to', 'ali_from', 'ali_to', 'description', 'taxon_id', 'enzyme_id',
#'   'source_file', 'hmm_coverage', 'target_coverage', 'hit_id'.
#'
#' @details
#' The function:
#' - reads each file line-by-line,
#' - skips comment lines beginning with '#',
#' - tokenizes fields by whitespace,
#' - constructs numeric columns with 'suppressWarnings(as.numeric(...))' to
#'   avoid hard failures on occasional non-numeric fields,
#' - computes 'hmm_coverage = (hmm_to - hmm_from + 1) / qlen' and
#'   'target_coverage = (ali_to - ali_from + 1) / tlen',
#' - attaches 'taxon_id' and 'enzyme_id' inferred from filename.
#'
#'
#' @examples
#' # files <- list.files("hmmsearch", pattern = "\\.(domtbl|tbl)$", full.names = TRUE)
#' # hits_df <- parse_hmmsearch_tbl(files)
#'
#' @export
parse_hmmsearch_tbl <- function(files) {
  
  # accumulator for per-file data.frames
  all_hits <- list()
  n <- 0L
  
  # iterate files
  for (f in files) {
    # infer metadata from filename, skip files that don't match convention
    meta <- tryCatch(parse_hmm_filename(f), error = function(e) NULL)
    if (is.null(meta)) next
    
    # read lines and remove comments/empty lines
    lines <- readLines(f, warn = FALSE)
    lines <- lines[!grepl("^#", lines) & nzchar(lines)]
    if (length(lines) == 0) next
    
    # parse each valid line to a small data.frame
    parsed <- lapply(lines, function(line) {
      # split tokens on whitespace
      toks <- strsplit(line, "\\s+", perl = TRUE)[[1]]
      # safety: skip malformed lines
      if (length(toks) < 18) return(NULL)
      
      n_toks <- length(toks)
      # description often starts at token 23 (if present)
      desc <- if (n_toks >= 23) paste(toks[23:n_toks], collapse = " ") else NA_character_
      
      # create a 1-row data.frame with parsed fields
      data.frame(
        target_name = toks[1],
        tlen = suppressWarnings(as.numeric(toks[3])),
        query_name = toks[4],
        qlen = suppressWarnings(as.numeric(toks[6])),
        full_evalue = suppressWarnings(as.numeric(toks[7])),
        full_score = suppressWarnings(as.numeric(toks[8])),
        dom_num = suppressWarnings(as.numeric(toks[10])),
        dom_of = suppressWarnings(as.numeric(toks[11])),
        i_evalue = suppressWarnings(as.numeric(toks[13])),
        dom_score = suppressWarnings(as.numeric(toks[14])),
        dom_bias = suppressWarnings(as.numeric(toks[15])),
        hmm_from = suppressWarnings(as.numeric(toks[16])),
        hmm_to = suppressWarnings(as.numeric(toks[17])),
        ali_from = suppressWarnings(as.numeric(toks[18])),
        ali_to = suppressWarnings(as.numeric(toks[19])),
        description = desc,
        stringsAsFactors = FALSE
      )
    })
    
    # drop NULLs (failed lines)
    parsed <- parsed[!vapply(parsed, is.null, logical(1))]
    if (length(parsed) == 0) next
    
    # combine rows for this file
    df <- do.call(rbind, parsed)
    
    # attach metadata derived from filename
    df$taxon_id <- meta$taxon_id
    df$enzyme_id <- meta$enzyme_id
    df$source_file <- basename(f)
    
    # compute coverage metrics
    df$hmm_coverage <- with(df, (hmm_to - hmm_from + 1) / qlen)
    df$target_coverage <- with(df, (ali_to - ali_from + 1) / tlen)
    
    # generate a unique hit identifier
    df$hit_id <- paste(df$taxon_id, df$enzyme_id, df$target_name, df$dom_num, sep = "|")
    
    # append
    n <- n + 1L
    all_hits[[n]] <- df
  }
  
  # combine all per-file data.frames, or return empty df if nothing parsed
  if (n == 0L) return(data.frame())
  res <- do.call(rbind, all_hits)
  rownames(res) <- NULL
  res
}
