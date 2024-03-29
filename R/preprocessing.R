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
  if (length(setdiff(samples, phenotype$table$Experiment)) > 0) {
    stop("The phenotypes of some samples were missing!")
  }
  new.object = data.table::copy(object)
  d.f = as_data_frame(phenotype$table, row.names = "Experiment")
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
  if (length(setdiff(samples, phenotype$table$Experiment)) > 0) {
    stop("The phenotypes of some samples were missing!")
  }
  new.object = data.table::copy(object)
  d.f = as_data_frame(phenotype$table, row.names = "Experiment")
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
    contaminant = grep("contaminant", names(new.object[[i]]), ignore.case = TRUE)
    reverse = grep("Reverse", names(new.object[[i]]), ignore.case = TRUE)
    only_by_site = grep("Only identified by site", names(new.object[[i]]), ignore.case = TRUE)
    if (length(contaminant) == 0 && length(reverse) == 0 && length(only_by_site) == 0) {
        warning("No column for checking false hits!")
    } else {
      if (length(contaminant) > 0) {
        contaminant = assert_length_1(contaminant)
        new.object[[i]] = new.object[[i]][
          is.na(new.object[[i]][[contaminant]]) | new.object[[i]][[contaminant]] != "+", ]
      }
      if (length(reverse) > 0) {
        reverse = assert_length_1(reverse)
        new.object[[i]] = new.object[[i]][
          is.na(new.object[[i]][[reverse]]) | new.object[[i]][[reverse]] != "+", ]
      }
      if (length(only_by_site) > 0) {
        only_by_site = assert_length_1(only_by_site)
        new.object[[i]] = new.object[[i]][
          is.na(new.object[[i]][[only_by_site]]) | new.object[[i]][[only_by_site]] != "+", ]
      }
    }
  }
  feature.counts[, "Remove false hits" := vapply(new.object, nrow, integer(1))]
  ## Get gene names
  get_genename_from_fasta_header = function(x) {
    ifelse(grepl("OS=Homo .+GN=", x),
      yes = strsplit(x, "GN=") %>%
        lapply(function(v) sub("\\b.*$", "", v[-1])) %>%
        vapply(function(v) paste(unique(v), collapse = ";"), character(1)),
      no = character(1)
      )
  }
  for (i in 1:length(new.object)) {
    if ("Gene names" %in% names(new.object[[i]])) {
      if ("Fasta headers" %in% names(new.object[[i]])) {
        genes = get_genename_from_fasta_header(new.object[[i]][["Fasta headers"]])
        new.object[[i]][, "Gene names" := ifelse(nzchar(genes), genes, new.object[[i]][["Gene names"]])]
      }
      new.object[[i]] = new.object[[i]][nzchar(new.object[[i]][["Gene names"]]), ]
    } else {
      if ("Fasta headers" %in% names(new.object[[i]])) {
        new.object[[i]][, "Gene names" := new.object[[i]][["Fasta headers"]]]
      }
    }
  }
  feature.counts[, "With gene names" := vapply(new.object, nrow, integer(1))]
  ## Check unique peptides
  for (i in 1:length(new.object)) {
    unique_peptides = grep("^Unique peptides", names(new.object[[i]]), ignore.case = TRUE)
    if (length(unique_peptides) < 1) {
      warning("No column for checking unique peptides!")
    } else {
      unique_peptides = assert_length_1(unique_peptides)
      new.object[[i]] = new.object[[i]][new.object[[i]][[unique_peptides]] >= min.unique.peptides, ]
    }
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
#' 
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
  if (length(inverse) != 1 && length(inverse) != length(object)) {
    stop("Invalid parameter 'inverse'!")
  }
  transform = function(x, opposite) {
    if (opposite) -log2(x) else log2(x)
  }
  new.object = data.table::copy(object)
  new.object[1:length(new.object)] = 
    mapply(transform, new.object, inverse, SIMPLIFY = FALSE)
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
#' @param method (Character) Choose the method(s) to use.
#' @return A new object of class ExpAssayFrame.
#' 
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
  if (length(method) != 1 && length(method) != length(object)) {
    stop("Invalid parameter 'method'!")
  }
  if (!all(method %in% c("center.median", "zscore"))) {
    stop("The input method for normalization was not supported!")
  }
  calculate = list(
    "center.median" = function(x) {
      x[is.na(x) | is.infinite(x)] = NA
      x - stats::median(x, na.rm = TRUE)
    }, 
    "zscore" = function(x) {
      x[is.na(x) | is.infinite(x)] = NA
      (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)
    }
  )
  normalize = function(x, method) {
    d.f = as.data.frame(t(x))
    d.f[, colnames(d.f)] = lapply(d.f, calculate[[method]])
    as.data.frame(t(d.f))
  }
  new.object = data.table::copy(object)
  new.object[1:length(new.object)] = 
    mapply(normalize, new.object, method, SIMPLIFY = FALSE)
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

