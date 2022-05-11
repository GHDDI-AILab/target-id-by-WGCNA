#' @importFrom data.table fread setnames
#' @importFrom utils str
#' 
NULL

#' Read the experiment information of mass spectrometry data.
#'
#' Load the data from the experimentalDesign.txt file in MaxQuantOutput and create an ExperimentInfo object.
#'
#' @param data.dir Directory containing the MaxQuantOutput data.
#' @return A list containing the experiment information.
#' @export
#' @examples
#' \dontrun{
#' object <- ReadExperimentalDesign('.')
#' }
#'
ReadExperimentalDesign = function(
  data.dir = '.'
) {
  ## Check the infile
  pattern = c(
    'experimentalDesign*.*',
    'MaxQuantOutput*/experimentalDesign*.*'
    )
  infile = Sys.glob(file.path(data.dir, pattern))
  if (length(infile) < 1) {
    stop('No experimentalDesign file was found!')
  } else if (length(infile) > 1) {
    warning('More than one experimentalDesign file!')
  }
  ## Load the table
  DT = data.table::fread(infile[[1]])
  if (! nrow(DT) || ! 'Experiment' %in% names(DT)) {
    stop('Empty or invalid experimentalDesign file!')
  }
  structure(
    list(experiments = unique(sort(DT$Experiment)), table = DT),
    filename = normalizePath(infile[[1]]),
    class = c('ExperimentInfo', 'ExperimentList', 'list')
    )
}

#' Read the protein-level mass spectrometry data.
#'
#' Load the data from the proteinGroups.txt file in MaxQuantOutput and create an ExperimentAssay object.
#'
#' @param data.dir Directory containing the MaxQuantOutput data.
#' @param column.prefix A character vector with the prefixes of the expression columns to extract.
#' @return A list of data.tables containing the expression data.
#' @export
#' @examples
#' \dontrun{
#' object <- ReadProteinGroups('.', col = 'Ratio H/L normalized')
#' }
#'
ReadProteinGroups = function(
  data.dir = '.',
  column.prefix = c('Ratio H/L normalized', 'LFQ intensity')
) {
  ## Check the infile
  pattern = c(
    'proteinGroups*.*',
    'MaxQuantOutput*/proteinGroups*.*'
    )
  infile = Sys.glob(file.path(data.dir, pattern))
  if (length(infile) < 1) {
    stop('No proteinGroups file was found!')
  } else if (length(infile) > 1) {
    warning('More than one proteinGroups file!')
  }
  ## Load the tables
  info = ReadExperimentalDesign(data.dir = dirname(infile[[1]]))
  DT = data.table::fread(infile[[1]])
  if (! nrow(DT) || ! ncol(DT)) {
    stop('Empty or invalid proteinGroups file!')
  }
  ## Find the columns that are not sample-specific
  sample_cols = lapply(info$experiments, function(i) grep(i, names(DT)))
  non_samples = setdiff(1:ncol(DT), unlist(sample_cols))
  ## Find the columns that are sample-specific with the given prefixes
  Assays = list()
  for (prefix in column.prefix) {
    sample_cols = lapply(paste(prefix, info$experiments),
                         function(i) which(i == names(DT)))
    if (length(unlist(sample_cols))) {
      Assays[[prefix]] = cbind(
        DT[, non_samples, with = FALSE],
        DT[, unlist(sample_cols), with = FALSE]
        )
      data.table::setnames(Assays[[prefix]],
        old = paste(prefix, info$experiments), new = info$experiments,
        skip_absent = TRUE
        )
    }
  }
  rm(DT)
  structure(
    Assays,
    experiments = info$experiments,
    filename = normalizePath(infile[[1]]),
    class = c('ExperimentAssay', 'list')
    )
}

#' @rdname print
#' @method print ExperimentList
#'
print.ExperimentInfo = function(x) {
  cat(sprintf('An object of class %s\n\n', class(x)[1]))
  utils::str(x)
  invisible(x)
}

#' @rdname print
#' @method print ExperimentTable
#'
print.ExperimentAssay = function(x) {
  cat(sprintf('An object of class %s\n\n', class(x)[1]))
  cat(sprintf('%d assay(s):\n', length(x)))
  for (i in 1:length(x)) {
    cat(sprintf('\t%d features across %d samples within assay %d.\n',
                nrow(x[[i]]), length(attr(x, 'experiments')), i))
  }
  invisible(x)
}

