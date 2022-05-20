#' Print the name and the value of a variable.
#' 
#' @param x The object to display.
#' @return The input variable.
#' @examples
#' \dontrun{
#' data.dir = '/path/to/dir'
#' Show(data.dir)
#' }
#' 
Show = function(x) {
  name = as.character(substitute(x))
  repr = capture.output(dput(x))
  cat(sprintf("%s = %s\n", name, repr))
  invisible(x)
}

#' Assert a length-1 vector.
#' 
#' Execute an error if the length of the input is 0 and 
#' execute a warning if the length of the input is bigger 
#' than 1. Return the first element of the input.
#' 
#' @param x The object for checking.
#' @return The first element of the input.
#' @examples
#' \dontrun{
#' data.dir = '.'
#' Show(data.dir)
#' assert_length_1(data.dir)
#' data.dir = c('.', '/path/to/dir')
#' Show(data.dir)
#' assert_length_1(data.dir)
#' }
#' 
assert_length_1 = function(x) {
  if (length(x) < 1) {
    name = as.character(substitute(x))
    stop(sprintf('No %s was found!', name))
  } else if (length(x) > 1) {
    name = as.character(substitute(x))
    warning(sprintf('More than one %s!', name))
  }
  return(x[1])  ## `[` won't change the type. 
}

#' Get the gene information from org.Hs.eg.db.
#' 
get_geneinfo = function(
  d.t, 
  remove_dup = TRUE
) {

#' A dict with previous and current gene names.
GENES_DICT = c(
  # previous / approved names
  "AARS"   = "AARS1", 
  "ATP5A1" = "ATP5F1A", 
  "ATP5B"  = "ATP5F1B", 
  "ATP5F1" = "ATP5PB", 
  "GARS"   = "GARS1", 
  "H1F0"   = "H1-0", 
  "WARS"   = "WARS1"
  )

  d.t[gene %in% names(GENES_DICT), gene := GENES_DICT[gene]]
  select(org.Hs.eg.db::org.Hs.eg.db, 
    keys = d.t[, gene], columns = c("ENSEMBL","GENENAME"), keytype = "SYMBOL"
    ) %>% 
  suppressMessages() %>% 
  data.table::as.data.table() %>% 
  data.table::setnames(., 
    c("ENSEMBL", "GENENAME", "SYMBOL"), 
    c("ensembl", "fullname", "gene")) %>% 
  {if (remove_dup) .[, head(.SD, 1), by = "gene"] else .} %>% 
  .[d.t, on = "gene"]
}

