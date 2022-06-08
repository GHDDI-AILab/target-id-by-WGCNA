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
    stop(sprintf('No %s was found!', as.character(substitute(x))))
  } else if (length(x) > 1) {
    warning(sprintf('More than one %s!', as.character(substitute(x))))
  }
  return(x[1])  ## `[` won't change the type. 
}

#' Get the gene information from org.Hs.eg.db.
#' 
#' @param d.t A data.table.
#' @param remove_dup (A length-1 logical) Remove duplicate genes or not.
#' @return A data.table.
#' @examples
#' \dontrun{
#' d.t = data.table::data.table(gene = "ATP5A1")
#' d.t = get_geneinfo(d.t)
#' }
#' 
get_geneinfo = function(
  d.t, 
  remove_dup = TRUE
) {

  ## A dict with previous and current gene names.
  GENES_DICT = c(
    # previous / approved names
    "AARS"   = "AARS1", 
    "ATP5A1" = "ATP5F1A", 
    "ATP5B"  = "ATP5F1B", 
    "ATP5F1" = "ATP5PB", 
    "ATP5O"  = "ATP5PO", 
    "CSDA"   = "YBX3", 
    "GARS"   = "GARS1", 
    "H1F0"   = "H1-0", 
    "WARS"   = "WARS1"
    )

  d.t[gene %in% names(GENES_DICT), gene := GENES_DICT[gene]]
  AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, 
    keys = d.t[["gene"]], columns = c("ENSEMBL","GENENAME"), keytype = "SYMBOL"
    ) %>% 
  suppressMessages() %>% 
  data.table::as.data.table() %>% 
  data.table::setnames(., 
    c("ENSEMBL", "GENENAME", "SYMBOL"), 
    c("ensembl", "fullname", "gene")) %>% 
  {if (remove_dup) .[, head(.SD, 1), by = "gene"] else .} %>% 
  .[d.t, on = "gene"]
}

#' Check missing values in mass spectrometry data.
#' 
#' Check if the values are missing in MS data 
#' (NA in labelled MS data, and 0 in label-free MS data).
#' @param x A numeric vector.
#' @return A logical vector.
#' @examples
#' \dontrun{
#' is_missing_in_ms(c(1, 0, NA, -Inf))
#' }
#' 
is_missing_in_ms = function(x) {
  is.na(x) | x == 0 | is.infinite(x)
}

