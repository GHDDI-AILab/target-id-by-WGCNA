#' Print the name and the value of a variable.
#' 
#' @param x The object to display.
#' @return The input variable.
#' @export
#' @examples
#' data.dir = "/path/to/dir"
#' Show(data.dir)
#' 
Show = function(x) {
  name = capture.output(dput(substitute(x)))
  repr = capture.output(dput(x))
  cat(sprintf("%s = %s\n", name, repr))
  invisible(x)
}

#' Assert a length-1 vector or list.
#' 
#' Execute an error if the length of the input is 0 and 
#' execute a warning if the length of the input is bigger 
#' than 1. Return a new object of the same class as the 
#' input, with the first element of the input.
#' 
#' @param x The object for checking.
#' @return An object with the first element of the input.
#' @export
#' @examples
#' data.dir = "."
#' Show(data.dir)
#' assert_length_1(data.dir)
#' data.dir = c(".", "/path/to/dir")
#' Show(data.dir)
#' assert_length_1(data.dir)
#' 
assert_length_1 = function(x) {
  if (length(x) < 1) {
    stop(sprintf('No %s was found!', capture.output(dput(substitute(x)))))
  } else if (length(x) > 1) {
    warning(sprintf('More than one %s!', capture.output(dput(substitute(x)))))
  }
  return(x[1])  ## `[` won't change the type. 
}

#' Get the gene information from org.Hs.eg.db.
#' 
#' @param d.t A data.table.
#' @param remove_dup (A length-1 logical) Remove duplicate genes or not.
#' @return A data.table.
#' @export
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
  d.t = data.table::as.data.table(d.t)
  if (all(c("ensembl", "fullname", "gene") %in% names(d.t))) {
    return(d.t[])
  }
  ## A dict with previous and current gene names.
  GENES_DICT = c(
    ## previous / approved names
    "AARS"      = "AARS1", 
    "ADSS"      = "ADSS2", 
    "ATP5A1"    = "ATP5F1A", 
    "ATP5B"     = "ATP5F1B", 
    "ATP5F1"    = "ATP5PB", 
    "ATP5H"     = "ATP5PD", 
    "ATP5O"     = "ATP5PO", 
    "C19orf10"  = "MYDGF", 
    "CSDA"      = "YBX3", 
    "GARS"      = "GARS1", 
    "H1F0"      = "H1-0", 
    "H1FV"      = "H1-0", 
    "H1F1"      = "H1-1", 
    "H1F2"      = "H1-2", 
    "H1F3"      = "H1-3", 
    "H1F4"      = "H1-4", 
    "H1F5"      = "H1-5", 
    "H1FX"      = "H1-10", 
    "HIST1H1A"  = "H1-1", 
    "HIST1H1C"  = "H1-2", 
    "HIST1H1D"  = "H1-3", 
    "HIST1H1E"  = "H1-4", 
    "HIST1H1B"  = "H1-5", 
    "HIST1H1T"  = "H1-6", 
    "H1FT"      = "H1-6", 
    "H1FNT"     = "H1-7", 
    "H1FOO"     = "H1-8", 
    "H1-9"      = "H1-9P", 
    "HIST1H2BA" = "H2BC1", 
    "HIST1H2BB" = "H2BC3", 
    "HIST1H2BC" = "H2BC4", 
    "HIST1H2BD" = "H2BC5", 
    "HIST1H2BE" = "H2BC6", 
    "HIST1H2BF" = "H2BC7", 
    "HIST1H2BG" = "H2BC8", 
    "HIST1H2BH" = "H2BC9", 
    "HIST1H2BI" = "H2BC10", 
    "HIST1H2BJ" = "H2BC11", 
    "HIST1H2BK" = "H2BC12", 
    "HIST1H2BL" = "H2BC13", 
    "HIST1H2BM" = "H2BC14", 
    "HIST1H2BN" = "H2BC15", 
    "HIST1H2BO" = "H2BC17", 
    "KARS"      = "KARS1", 
    "PKM2"      = "PKM", 
    "PTRF"      = "CAVIN1", 
    "TOMM70A"   = "TOMM70", 
    "VARS"      = "VARS1", 
    "WARS"      = "WARS1"
    )
  ## Check the column `gene`
  if (! "gene" %in% names(d.t)) {
    col = grep("^gene", names(d.t), ignore.case = TRUE, value = TRUE) %>% 
      assert_length_1()
    setnames(d.t, col, "gene")
  } 
  d.t[gene %in% names(GENES_DICT), "gene" := GENES_DICT[gene]]
  AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, 
    keys = d.t[["gene"]], columns = c("ENSEMBL","GENENAME"), keytype = "SYMBOL"
    ) %>% 
  suppressMessages() %>% 
  data.table::as.data.table() %>% 
  data.table::setnames(
    c("ENSEMBL", "GENENAME", "SYMBOL"), 
    c("ensembl", "fullname", "gene")) %>% 
  {if (remove_dup) .[, head(.SD, 1), by = "gene"] else .} %>% 
  .[d.t, on = "gene"]
}

#' Check missing values in mass spectrometry data.
#' 
#' Check if the values are missing in MS data 
#' (NA in labelled MS data, and 0 in label-free MS data).
#' 
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

#' Coerce to a data.frame.
#' 
#' Coerce to a data.frame, and move a column 
#' to be row names if specifying the column.
#' 
#' @param x An object.
#' @param row.names A length-1 character specifying the column to be row names, 
#'   or a call, e.g. `x$ID` or `x[["ID"]]`, specifying the column. Else, this param 
#'   will be passed to as.data.frame().
#' @param ... additional arguments to be passed to as.data.frame().
#' @return A data.frame object.
#' @export
#' @examples
#' a = data.table::data.table(name = LETTERS[1:10], value = runif(10))
#' as_data_frame(a, a$name)
#' 
as_data_frame = function(
  x, 
  row.names = NULL, 
  ...
) {
  if (is.character(substitute(row.names)) && length(row.names) == 1 && !is.na(row.names)) {
    d_f = as.data.frame(x, ...)
    rownames(d_f) = d_f[[row.names]]
    d_f[[row.names]] = NULL
  } else if (is.call(substitute(row.names)) && length(substitute(row.names)) == 3) {
    d_f = as.data.frame(x, ...)
    name = as.character(substitute(row.names)[[3]])
    rownames(d_f) = d_f[[name]]
    d_f[[name]] = NULL
  } else {
    d_f = as.data.frame(x, row.names = row.names, ...)
  }
  return(d_f)
}

