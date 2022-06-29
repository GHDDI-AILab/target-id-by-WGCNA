#' @include util.R
#' @include generics.R
#' 
NULL

#' Tidy an ExpAssayTable object.
#' 
#' @param object An object of class ExpAssayTable.
#' @return A new object of class ExpAssayTable.
#' 
#' @rdname Tidy
#' @method Tidy ExpAssayTable
#' @export
#' @examples
#' \dontrun{
#' new.Assay = Tidy(old.Assay)
#' }
#' 
Tidy.ExpAssayTable = function(
  object
) {
  new.object = data.table::copy(object)
  samples = attr(new.object, "experiments")
  for (i in 1:length(new.object)) {
    ## Solve the "integer64" problem
    new.object[[i]][, (samples) := 
      lapply(new.object[[i]][, samples, with = FALSE], as.double)][]
    ## Convert NaN or 0 to NA
    convert = function(x) ifelse(is_missing_in_ms(x), NA, x)
    new.object[[i]][, (samples) := 
      lapply(new.object[[i]][, samples, with = FALSE], convert)][]
  }
  return(new.object)
}

#' @rdname Tidy
#' @method Tidy default
#' @export
Tidy.default = function(object, ...) {
  if (inherits(object, "ExpAssayTable")) {
    Tidy.ExpAssayTable(object, ...)
  } else {
    stop("This method is associated with class ExpAssayTable.")
  }
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
  new.samples = intersect(samples, attr(object, "experiments"))
  if (length(new.samples) < 1) {
    stop("Invalid input samples!")
  } else if (length(new.samples) != length(samples)) {
    warning("Some of the given sample names were invalid!")
  }
  new.object = data.table::copy(object)
  for (i in 1:length(new.object)) {
    non_samples = setdiff(names(new.object[[i]]), attr(new.object, "experiments"))
    new.object[[i]] = new.object[[i]][, c(non_samples, new.samples), with = FALSE]
  }
  attr(new.object, "experiments") = new.samples
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
  new.samples = intersect(samples, attr(object, "experiments"))
  if (length(new.samples) < 1) {
    stop("Invalid input samples!")
  } else if (length(new.samples) != length(samples)) {
    warning("Some of the given sample names were invalid!")
  }
  new.object = data.table::copy(object)
  for (i in 1:length(new.object)) {
    new.object[[i]] = new.object[[i]][new.samples, ]
  }
  attr(new.object, "experiments") = new.samples
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

#' Add a phenotype table of the samples in the assay.
#' 
#' @param object An object of class ExpAssayTable.
#' @param phenotype An object of class ExperimentInfo.
#' @return A new object of class ExpAssayTable.
#' 
#' @rdname AddPhenotype
#' @method AddPhenotype ExpAssayTable
#' @export
#' 
AddPhenotype.ExpAssayTable = function(
  object, 
  phenotype
) {
  ATTR_PHENO = "phenotype"
  samples = attr(object, "experiments")
  if (length(setdiff(samples, phenotype$experiments)) > 0) {
    stop("The phenotypes of some samples were missing!")
  }
  d.f = as.data.frame(phenotype$table)
  rownames(d.f) = d.f$Experiment
  d.f$Experiment = NULL
  new.object = data.table::copy(object)
  attr(new.object, ATTR_PHENO) = d.f
  return(new.object)
}

#' Add a phenotype table of the samples in the assay.
#' 
#' @param object An object of class ExpAssayFrame.
#' @param phenotype An object of class ExperimentInfo.
#' @return A new object of class ExpAssayFrame.
#' 
#' @rdname AddPhenotype
#' @method AddPhenotype ExpAssayFrame
#' @export
#' 
AddPhenotype.ExpAssayFrame = function(
  object, 
  phenotype
) {
  ATTR_PHENO = "phenotype"
  samples = attr(object, "experiments")
  if (length(setdiff(samples, phenotype$experiments)) > 0) {
    stop("The phenotypes of some samples were missing!")
  }
  d.f = as.data.frame(phenotype$table)
  rownames(d.f) = d.f$Experiment
  d.f$Experiment = NULL
  new.object = data.table::copy(object)
  attr(new.object, ATTR_PHENO) = d.f
  return(new.object)
}

#' @rdname AddPhenotype
#' @method AddPhenotype default
#' @export
AddPhenotype.default = function(object, ...) {
  if (inherits(object, "ExpAssayTable")) {
    AddPhenotype.ExpAssayTable(object, ...)
  } else if (inherits(object, "ExpAssayFrame")) {
    AddPhenotype.ExpAssayFrame(object, ...)
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
  min.fraction = 0.5
) {
  ## Create a new object
  new.object = data.table::copy(object)
  feature.counts = data.table::data.table()
  feature.counts[, "Assay" := names(new.object)]
  feature.counts[, "Raw data" := vapply(new.object, nrow, integer(1))]
  ## Remove false hits
  ## column name: `Potential contaminant` or `Contaminant`
  ## column name: `Reverse`
  ## column name: `Only identified by site`
  for (i in 1:length(new.object)) {
    d.t = new.object[[i]]
    contaminant  = grep("contaminant", names(d.t), 
                        ignore.case = TRUE, value = TRUE) %>% assert_length_1()
    reverse      = grep("Reverse", names(d.t), 
                        ignore.case = TRUE, value = TRUE) %>% assert_length_1()
    only_by_site = grep("Only identified by site", names(d.t), 
                        ignore.case = TRUE, value = TRUE) %>% assert_length_1()
    new.object[[i]] = d.t[
      (is.na(d.t[[contaminant]]) | d.t[[contaminant]] != "+") & 
      (is.na(d.t[[reverse]]) | d.t[[reverse]] != "+") & 
      (is.na(d.t[[only_by_site]]) | d.t[[only_by_site]] != "+"), 
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
    non.NA.fraction = rowSums(!is_missing_in_ms(mat)) / nSamples
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
#' Assay = ReadProteinGroups(".")
#' assay = Reshape(Assay)
#' }
#' 
Reshape.ProteinGroups = function(
  object
) {
  new.object = Tidy(object)
  samples = attr(new.object, "experiments")
  for (i in 1:length(new.object)) {
    genenames = new.object[[i]][["Gene names"]]
    new.object[[i]] = as.data.frame(t(new.object[[i]][, samples, with = FALSE]))
    colnames(new.object[[i]]) = genenames
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

#' Perform log2 transformation of an expression profile.
#' 
#' @param object An object of class ExpAssayFrame.
#' @param inverse (logical) Compute log2 (FALSE) or -log2 (TRUE).
#' @return A new object of class ExpAssayFrame.
#' @rdname LogTransform
#' @method LogTransform ExpAssayFrame
#' @export
#' @examples
#' \dontrun{
#' new.assay = LogTransform(old.assay)
#' }
#' 
LogTransform.ExpAssayFrame = function(
  object, 
  inverse = FALSE
) {
  new.object = data.table::copy(object)
  for (i in 1:length(new.object)) {
    new.object[[i]] = 
      if (inverse) -log2(new.object[[i]]) else log2(new.object[[i]])
  }
  return(new.object)
}

#' @rdname LogTransform
#' @method LogTransform default
#' @export
#' 
LogTransform.default = function(object, ...) {
  if (inherits(object, "ExpAssayFrame")) {
    LogTransform.ExpAssayFrame(object, ...)
  } else {
    stop("This method is associated with class ExpAssayFrame.")
  }
}

#' Normalize a log-transformed expression profile.
#' 
#' @param object An object of class ExpAssayFrame.
#' @param method Choose the method to use.
#' @return A new object of class ExpAssayFrame.
#' @rdname Normalize
#' @method Normalize ExpAssayFrame
#' @export
#' @examples
#' \dontrun{
#' new.assay = Normalize(LogTransform(old.assay), method = "center.median")
#' }
#' 
Normalize.ExpAssayFrame = function(
  object, 
  method = "center.median"
) {
  if (assert_length_1(method) == "center.median") {
    normalize = function(x) {
      x[is_missing_in_ms(x)] = NA
      x - stats::median(x, na.rm = TRUE)
    }
  } else {
    stop("The input method for normalization was not supported!")
  }
  new.object = data.table::copy(object)
  for (i in 1:length(new.object)) {
    d.f = as.data.frame(t(new.object[[i]]))
    d.f[, colnames(d.f)] = lapply(d.f, normalize)
    new.object[[i]] = as.data.frame(t(d.f))
  }
  return(new.object)
}

#' @rdname Normalize
#' @method Normalize default
#' @export
#' 
Normalize.default = function(object, ...) {
  if (inherits(object, "ExpAssayFrame")) {
    Normalize.ExpAssayFrame(object, ...)
  } else {
    stop("This method is associated with class ExpAssayFrame.")
  }
}

