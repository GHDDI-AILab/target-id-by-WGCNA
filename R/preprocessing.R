#' @include util.R
#' @include generics.R
#' 
NULL

#' Read the experiment information of mass spectrometry data.
#'
#' Load the data from the experimentalDesign.txt file in MaxQuantOutput/ and create an ExperimentInfo object.
#'
#' @param data.dir (A length-1 character) Directory containing the MaxQuantOutput data.
#' @return (An ExperimentInfo object) A list containing the experiment information.
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
  structure(
    list(experiments = unique(sort(DT$Experiment)), table = DT),
    filename = normalizePath(infile),
    class = c('ExperimentInfo', 'ExperimentList', 'list')
    )
}

#' Read the protein-level mass spectrometry data.
#'
#' Load the data from the proteinGroups.txt file in MaxQuantOutput/ and create a ProteinGroups object.
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
  sample_cols = lapply(info$experiments, function(i) grep(i, names(DT)))
  non_samples = setdiff(1:ncol(DT), unlist(sample_cols))
  ## Find the columns that are sample-specific with the given prefixes
  Assays = list()
  for (prefix in column.prefix) {
    ## There may be no column with the prefix, so we need to check:
    sample_cols = lapply(paste(prefix, info$experiments),
                         function(i) which(i == names(DT))) %>% unlist
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
  rm(DT)
  structure(
    Assays,
    experiments = info$experiments,
    filename = normalizePath(infile),
    class = c('ProteinGroups', 'ExperimentAssay', 'list')
    )
}

#' @rdname print
#' @method print ExperimentList
#' @export
#'
print.ExperimentList = function(x) {
  cat(sprintf('An object of class %s\n\n', class(x)[1]))
  utils::str(x)
  invisible(x)
}

#' @rdname print
#' @method print ExperimentAssay
#' @export
#'
print.ExperimentAssay = function(x) {
  samples = attr(x, 'experiments')
  cat(sprintf('An object of class %s\n\n', class(x)[1]))
  cat(sprintf('%d Experiment(s): ', length(samples)))
  if (length(samples) > 3) {
    cat(sprintf('"%s", "%s", "%s", ...\n', 
                samples[1], samples[2], samples[3]))
  } else {
    cat(sprintf('"%s"\n', paste(samples, collapse = '", "')))
  }
  cat(sprintf('%d Assay(s): "%s"\n', 
              length(x), paste(names(x), collapse = '", "')))
  for (i in 1:length(x)) {
    cat(sprintf('\t%d features across %d samples within assay %d.\n',
                nrow(x[[i]]), length(attr(x, 'experiments')), i))
  }
  invisible(x)
}

#' Subsetting an ExperimentAssay object.
#' 
#' @param object An object of class ExperimentAssay.
#' @param samples (A character vector) The samples for subsetting.
#' @return A new object of class ExperimentAssay.
#' 
#' @rdname Subset
#' @method Subset ExperimentAssay
#' @examples
#' \dontrun{
#' new.Assay = Subset(old.Assay, samples = c("sample_1", "sample_2", "sample_3"))
#' }
#' 
Subset.ExperimentAssay = function(object, samples) {
  if (!is.character(samples)) {
    stop("Invalid input samples!")
  }
  new.object = data.table::copy(object)
  experiments = attr(new.object, "experiments")
  for (i in 1:length(new.object)) {
    non_samples = setdiff(names(new.object[[i]]), experiments)
    new.object[[i]] = new.object[[i]][, c(non_samples, samples), with = FALSE]
  }
  attr(new.object, "experiments") = samples
  return(new.object)
}

#' @rdname Subset
#' @method Subset default
#' @export
Subset.default = function(object, ...) {
  if (inherits(object, "ExperimentAssay")) {
    Subset.ExperimentAssay(object, ...)
  } else {
    stop("This method is associated with class ExperimentAssay.")
  }
}

#' Quality control of MS-based proteomics data.
#' 
#' @param object An object of class ProteinGroups.
#' @param min.unique.peptides A threshold for the number of unique peptides.
#' @param min.fraction A minimum fraction of non-missing samples for a gene 
#'    to be considered good.
#' @return A new object of class ProteinGroups.
#' 
#' @rdname QC
#' @method QC ProteinGroups
#' @examples
#' \dontrun{
#' new.Assay = QC(old.Assay)
#' }
#' 
QC.ProteinGroups = function(
  object, 
  min.unique.peptides = 2, 
  min.fraction = 0.5
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
        lapply(., function(v) sub(" .*$", "", v[-1])) %>%
        vapply(., function(v) paste(unique(v), collapse = ";"), character(1)),
      no = character(1)
      )
  }
  for (i in 1:length(new.object)) {
    genes = new.object[[i]][, get_genename_from_fasta_header(`Fasta headers`)]
    new.object[[i]][, `Gene names` := ifelse(nzchar(genes), genes, `Gene names`)]
    new.object[[i]] = new.object[[i]][nzchar(`Gene names`), ]
  }
  feature.counts[, "With gene names" := vapply(new.object, nrow, integer(1))]
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
    non.NA.fraction = rowSums(!is.na(mat)) / nSamples
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
    stop("This method is associated with class ExperimentAssay.")
  }
}

