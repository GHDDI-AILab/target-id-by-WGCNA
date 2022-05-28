#' @include util.R
#' @include generics.R
#' 
NULL

#' Read the experiment information of mass spectrometry data.
#'
#' Load the data from the experimentalDesign.txt file in 
#' MaxQuantOutput/ and create an ExperimentInfo object.
#'
#' @param data.dir (length-1 character vector) 
#'   Directory containing the MaxQuantOutput data.
#' @return (An ExperimentInfo object) 
#'   A list containing the experiment information.
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
  data.dir = assert_length_1(data.dir)
  experimentalDesign.txt = Sys.glob(file.path(data.dir, pattern))
  infile = assert_length_1(experimentalDesign.txt)
  ## Load the table
  DT = data.table::fread(infile)
  if (! nrow(DT) || ! 'Experiment' %in% names(DT)) {
    stop('Empty or invalid experimentalDesign file!')
  }
  ## Create an object
  structure(
    list(experiments = unique(sort(DT$Experiment)), table = DT),
    filename = normalizePath(infile),
    class = c('ExperimentInfo', 'ExperimentList', 'list')
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
  data.dir = assert_length_1(data.dir)
  proteinGroups.txt = Sys.glob(file.path(data.dir, pattern))
  infile = assert_length_1(proteinGroups.txt)
  ## Load the tables
  info = ReadExperimentalDesign(data.dir = dirname(infile))
  DT = data.table::fread(infile)
  if (! nrow(DT) || ! ncol(DT)) {
    stop('Empty or invalid proteinGroups file!')
  }
  ## Find the columns that are not sample-specific
  sample_cols = lapply(info$experiments, 
                       function(i) grep(i, names(DT))) %>% unlist()
  non_samples = setdiff(1:ncol(DT), sample_cols)
  ## Find the columns that are sample-specific with the given prefixes
  Assays = list()
  for (prefix in column.prefix) {
    ## There may be no column with the prefix, so we need to check:
    sample_cols = lapply(paste(prefix, info$experiments),
                         function(i) which(i == names(DT))) %>% unlist()
    if (length(sample_cols)) {
      Assays[[prefix]] = cbind(
        DT[, non_samples, with = FALSE],
        DT[, sample_cols, with = FALSE]
        )
      data.table::setnames(Assays[[prefix]],
        old = paste(prefix, info$experiments), new = info$experiments,
        skip_absent = TRUE
        )
    }
  }
  ## Remove the original data table from memory
  rm(DT)
  ## Create an object
  structure(
    Assays,
    experiments = info$experiments,
    filename = normalizePath(infile),
    class = c('ProteinGroups', 'ExpAssayTable', 'ExperimentList', 'list')
    )
}

#' Subsetting an ExpAssayTable object.
#' 
#' @param object An object of class ExpAssayTable.
#' @param samples (character) The samples for subsetting.
#' @return A new object of class ExpAssayTable.
#' 
#' @rdname Subset
#' @method Subset ExpAssayTable
#' @export
#' @examples
#' \dontrun{
#' new.Assay = Subset(old.Assay, samples = c("sample_1", "sample_2", "sample_3"))
#' }
#' 
Subset.ExpAssayTable = function(
  object, 
  samples
) {
  if (!is.character(samples)) {
    stop("Invalid input samples!")
  }
  new.object = data.table::copy(object)
  for (i in 1:length(new.object)) {
    non_samples = setdiff(names(new.object[[i]]), attr(new.object, "experiments"))
    new.object[[i]] = new.object[[i]][, c(non_samples, samples), with = FALSE]
  }
  attr(new.object, "experiments") = samples
  return(new.object)
}

#' Subsetting an ExpAssayFrame object.
#' 
#' @param object An object of class ExpAssayFrame.
#' @param samples (character) The samples for subsetting.
#' @return A new object of class ExpAssayFrame.
#' 
#' @rdname Subset
#' @method Subset ExpAssayFrame
#' @export
#' @examples
#' \dontrun{
#' new.assay = Subset(old.assay, samples = c("sample_1", "sample_2", "sample_3"))
#' }
#' 
Subset.ExpAssayFrame = function(
  object, 
  samples
) {
  if (!is.character(samples)) {
    stop("Invalid input samples!")
  }
  new.object = data.table::copy(object)
  for (i in 1:length(new.object)) {
    new.object[[i]] = new.object[[i]][samples, ]
  }
  attr(new.object, "experiments") = samples
  return(new.object)
}

#' @rdname Subset
#' @method Subset default
#' @export
Subset.default = function(object, ...) {
  if (inherits(object, "ExpAssayTable")) {
    Subset.ExpAssayTable(object, ...)
  } else if (inherits(object, "ExpAssayFrame")) {
    Subset.ExpAssayFrame(object, ...)
  } else {
    stop("This method is associated with class ExpAssayTable and ExpAssayFrame.")
  }
}

#' Perform quality control of MS-based proteomics data.
#' 
#' @param object An object of class ProteinGroups.
#' @param min.unique.peptides (numeric) A threshold for the number 
#'   of unique peptides.
#' @param min.fraction (numeric) A minimum fraction of non-missing 
#'   samples for a gene to be considered good.
#' @return A new object of class ProteinGroups.
#' 
#' @rdname QC
#' @method QC ProteinGroups
#' @export
#' @examples
#' \dontrun{
#' new.Assay = QC(old.Assay)
#' }
#' 
QC.ProteinGroups = function(
  object, 
  min.unique.peptides = 2, 
  min.fraction = 0.2
) {
  ## Create a new object
  new.object = data.table::copy(object)
  feature.counts = data.table::data.table()
  feature.counts[, "Assay" := names(new.object)]
  feature.counts[, "Raw data" := vapply(new.object, nrow, integer(1))]
  ## Remove false hits
  for (i in 1:length(new.object)) {
    new.object[[i]] = new.object[[i]][
      `Potential contaminant` != "+" & 
      `Reverse` != "+" & 
      `Only identified by site` != "+", 
      ]
  }
  feature.counts[, "Remove false hits" := vapply(new.object, nrow, integer(1))]
  ## Get gene names
  get_genename_from_fasta_header = function(x) {
    ifelse(grepl("OS=Homo .+GN=", x),
      yes = strsplit(x, "GN=") %>%
        lapply(., function(v) sub("\\b.*$", "", v[-1])) %>%
        vapply(., function(v) paste(unique(v), collapse = ";"), character(1)),
      no = character(1)
      )
  }
  if ("Gene names" %in% names(new.object[[1]])) {
    for (i in 1:length(new.object)) {
      genes = new.object[[i]][, get_genename_from_fasta_header(`Fasta headers`)]
      new.object[[i]][, "Gene names" := ifelse(nzchar(genes), genes, `Gene names`)]
      new.object[[i]] = new.object[[i]][nzchar(`Gene names`), ]
    }
    feature.counts[, "With gene names" := vapply(new.object, nrow, integer(1))]
  } else {
    for (i in 1:length(new.object)) {
      new.object[[i]][, "Gene names" := `Fasta headers`]
    }
  }
  ## Check unique peptides
  for (i in 1:length(new.object)) {
    new.object[[i]] = new.object[[i]][`Unique peptides` >= min.unique.peptides, ]
  }
  feature.counts[, paste("Unique peptides >=", min.unique.peptides) := 
		 vapply(new.object, nrow, integer(1))]
  ## Check non-missing fraction
  experiments = attr(new.object, "experiments")
  for (i in 1:length(new.object)) {
    nSamples = length(experiments)
    mat = as.matrix(new.object[[i]][, experiments, with = FALSE])
    non.NA.fraction = rowSums(!is_missing_value(mat)) / nSamples
    new.object[[i]] = new.object[[i]][non.NA.fraction >= min.fraction, ]
  }
  feature.counts[, paste("goodGenes, min.fraction >=", min.fraction) := 
		 vapply(new.object, nrow, integer(1))]
  attr(new.object, "QC") = feature.counts[]
  return(new.object)
}

#' @rdname QC
#' @method QC default
#' @export
QC.default = function(object, ...) {
  if (inherits(object, "ProteinGroups")) {
    QC.ProteinGroups(object, ...)
  } else {
    stop("This method is associated with class ExpAssayTable.")
  }
}

#' Reshape tables of an expression assay.
#' 
#' Extract only the expression levels and reshape the rows and 
#' columns for further analysis.
#' 
#' @param object An object of class ProteinGroups.
#' @return An object of class ExpAssayFrame, whose rows correspond 
#'   to samples and columns to genes.
#' 
#' @rdname Reshape
#' @method Reshape ProteinGroups
#' @export
#' @examples
#' \dontrun{
#' Assay = ReadProteinGroups('.')
#' assay = Reshape(Assay)
#' }
#' 
Reshape.ProteinGroups = function(
  object
) {
  new.object = data.table::copy(object)
  for (i in 1:length(new.object)) {
    ## TODO: the "integer64" problem
    new.object[[i]] = object[[i]][, attr(object, "experiments"), with = FALSE
				  ] %>% t() %>% as.data.frame()
    colnames(new.object[[i]]) = object[[i]][["Gene names"]]
  }
  class(new.object) = c("ExpAssayFrame", "ExperimentList", "list")
  return(new.object)
}

#' @rdname Reshape
#' @method Reshape default
#' @export
#' 
Reshape.default = function(object, ...) {
  if (inherits(object, "ProteinGroups")) {
    Reshape.ProteinGroups(object, ...)
  } else {
    stop("This method is associated with class ProteinGroups.")
  }
}

#' Perform log normalization of an expression profile.
#' 
#' @param object An object of class ExpAssayFrame.
#' @param inverse (logical) Compute log2 (FALSE) or -log2 (TRUE).
#' @return A new object of class ExpAssayFrame.
#' @rdname LogNorm
#' @method LogNorm ExpAssayFrame
#' @export
#' @examples
#' \dontrun{
#' new.assay = LogNorm(old.assay)
#' }
#' 
LogNorm.ExpAssayFrame = function(
  object, 
  inverse = TRUE
) {
  new.object = data.table::copy(object)
  for (i in 1:length(new.object)) {
    new.object[[i]] = 
      if (inverse) -log2(new.object[[i]]) else log2(new.object[[i]])
  }
  return(new.object)
}

#' @rdname LogNorm
#' @method LogNorm default
#' @export
#' 
LogNorm.default = function(object, ...) {
  if (inherits(object, "ExpAssayFrame")) {
    LogNorm.ExpAssayFrame(object, ...)
  } else {
    stop("This method is associated with class ExpAssayFrame.")
  }
}

