#' @include util.R
#' 
NULL

#' Read the experiment information of mass spectrometry data.
#'
#' Load the data from the experimentalDesign.txt or summary.txt file in 
#' MaxQuantOutput/ and create an ExperimentInfo object.
#'
#' @param data.dir (length-1 character vector) 
#'   Directory containing the MaxQuantOutput data.
#' @return (An ExperimentInfo object) 
#'   A list containing the experiment information.
#' @export
#' @examples
#' \dontrun{
#' object = ReadExperimentalDesign(".")
#' }
#'
ReadExperimentalDesign = function(
  data.dir = "."
) {
  ## Check the infile
  pattern = c(
    "MaxQuantOutput*/experimentalDesign*.*", 
    "experimentalDesign*.*", 
    "MaxQuantOutput*/summary*.*", 
    "summary*.*"
    )
  data.dir = assert_length_1(data.dir)
  experimentalDesign.txt = Sys.glob(file.path(data.dir, pattern))
  infile = assert_length_1(experimentalDesign.txt)
  ## Load the table
  DT = data.table::fread(infile)
  if (nrow(DT) < 1 || ! "Experiment" %in% names(DT)) {
    stop("Empty or invalid experimentalDesign file!")
  }
  ## Create an object
  DT[, Experiment := as.character(Experiment)]
  samples = setdiff(unique(DT$Experiment), "")
  structure(
    list(experiments = samples, 
         table = DT), 
    filename = normalizePath(infile),
    class = c("ExperimentInfo", "ExperimentList", "list")
    )
}

#' Read the protein-level mass spectrometry data.
#'
#' Load the data from the proteinGroups.txt file in 
#' MaxQuantOutput/ and create a ProteinGroups object.
#'
#' @param data.dir (length-1 character vector) 
#'   Directory containing the MaxQuantOutput data.
#' @param column.prefix A character vector with 
#'   the prefixes of the expression columns to extract.
#' @return (An ProteinGroups object) 
#'   A list of data.tables containing the expression data.
#' @export
#' @examples
#' \dontrun{
#' object = ReadProteinGroups(".", col = "Ratio H/L normalized")
#' }
#'
ReadProteinGroups = function(
  data.dir = ".",
  column.prefix = c("Ratio H/L normalized", "LFQ intensity", "Intensity")
) {
  ## Check the infile
  pattern = c(
    "MaxQuantOutput*/proteinGroups*.*", 
    "proteinGroups*.*"
    )
  data.dir = assert_length_1(data.dir)
  proteinGroups.txt = Sys.glob(file.path(data.dir, pattern))
  infile = assert_length_1(proteinGroups.txt)
  ## Load experimentalDesign.txt and proteinGroups.txt
  info = ReadExperimentalDesign(data.dir = dirname(infile))
  DT = data.table::fread(infile)
  if (nrow(DT) < 1 || ncol(DT) < 1) {
    stop("Empty or invalid proteinGroups file!")
  }
  ## Find the columns that are not sample-specific
  sample_cols = lapply(info$experiments, function(i) grep(i, names(DT))) %>% unlist()
  non_samples = setdiff(1:ncol(DT), sample_cols)
  ## Find the columns that are sample-specific with the given prefixes
  Assays = list()
  for (prefix in column.prefix) {
    ## There may be no column with the prefix, so we need to check:
    col_names = paste(prefix, info$experiments, sep = " ")
    sample_cols = lapply(col_names, function(i) which(i == names(DT))) %>% unlist()
    if (length(sample_cols) < 1) {
      col_names = paste(prefix, info$experiments, sep = "_")
      sample_cols = lapply(col_names, function(i) which(i == names(DT))) %>% unlist()
    }
    if (length(sample_cols) > 0) {
      Assays[[prefix]] = cbind(
        DT[, non_samples, with = FALSE],
        DT[, sample_cols, with = FALSE]
        )
      data.table::setnames(Assays[[prefix]],
        old = col_names, 
        new = as.character(info$experiments),
        skip_absent = TRUE
        )
    }
  }
  ## Check the extraction
  if (length(Assays) < 1) {
    stop("Cannot find the columns with the given column.prefix!")
  }
  ## Remove the original data table from memory
  rm(DT)
  ## Create an object
  structure(
    Assays, 
    filename = normalizePath(infile), 
    experiments = as.character(info$experiments), 
    phenotype = data.frame(), 
    QC = data.table::data.table(), 
    powerEstimate = list(), 
    network = list(), 
    connectivity = list(), 
    class = c("ProteinGroups", "ExpAssayTable", "ExperimentList", "list")
    )
}

#' Read the phenotype information of the samples.
#' 
#' @param file A length-1 character specifying the file path.
#' @return An ExperimentInfo object.
#' @export
#' 
ReadPhenotypeTable = function(
  file
) {
  COLS = c("Experiment", "experiment", "Sample", "sample")
  file = assert_length_1(file)
  ## Load the table
  DT = data.table::fread(file)
  if (nrow(DT) < 1 || ! any(COLS %in% names(DT))) {
    stop("Empty or invalid phenotype file!")
  }
  ## Create an object
  sample_column = intersect(COLS, names(DT))
  sample_column = assert_length_1(sample_column)
  setnames(DT, sample_column, "Experiment")
  DT[, Experiment := as.character(Experiment)]
  samples = setdiff(unique(DT$Experiment), "")
  structure(
    list(experiments = samples, 
         table = DT[Experiment %in% samples, head(.SD, 1), by = "Experiment"]), 
    filename = normalizePath(file),
    class = c("ExperimentInfo", "ExperimentList", "list")
    )
}

