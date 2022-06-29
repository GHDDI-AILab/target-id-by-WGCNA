#' @include util.R
#' @include generics.R
#' 
NULL

#' Compute module-trait correlation and significance.
#' 
#' @param object A CorrelationNetwork object.
#' @param samples A character vector specifying the rows for analysis.
#' @param traits A character vector specifying the columns of traits for analysis.
#' @param prefix (A length-1 character) A prefix representing the disease.
#' @return A new CorrelationNetwork object.
#' @export
#' 
ModuleSignificance.CorrelationNetwork = function(
  object, 
  samples, 
  traits, 
  prefix
) {
  ATTR_MS = "module-trait"
  ATTR_NET = "network"
  ATTR_PHENO = "phenotype"
  if (length(object) != length(attr(object, ATTR_NET))) {
    stop("Invalid CorrelationNetwork object in the input!")
  }
  samples = attr(object, "experiments")
  if (length(setdiff(samples, rownames(attr(object, ATTR_PHENO)))) > 0) {
    stop("The phenotypes of some samples were missing!")
  }
  if (missing(prefix)) {
    prefix = ""
  } else {
    prefix = paste0(assert_length_1(prefix), ".")
  }
  if (missing(traits) || length(traits) < 1) {
    traits = colnames(attr(object, ATTR_PHENO))
  }
  new.object = data.table::copy(object)
  attr(new.object, ATTR_MS) = list()
  for (i in 1:length(new.object)) {
    if (missing(samples)) {
      MEs = attr(object, ATTR_NET)[[i]][["moduleEigengenes"]]
      datTraits = attr(object, ATTR_PHENO)[rownames(MEs), traits]
    } else if (length(samples) > 1) {
      MEs = attr(new.object, ATTR_NET)[[i]][["moduleEigengenes"]][samples, ]
      datTraits = attr(object, ATTR_PHENO)[samples, traits]
    } else {
      stop("The number of samples should be more than one!")
    }
    module_trait_cor = WGCNA::cor(MEs, datTraits, use = "p")
    module_trait_pval = WGCNA::corPvalueStudent(module_trait_cor, nSamples = nrow(MEs))
    colnames(module_trait_cor) = paste0(prefix, colnames(module_trait_cor))
    colnames(module_trait_pval) = paste0("p.", prefix, colnames(module_trait_pval))
    attr(new.object, ATTR_MS)[[i]] = list(
      cor = as.data.frame(module_trait_cor), 
      pval = as.data.frame(module_trait_pval)
      )
  }
  names(attr(new.object, ATTR_MS)) = names(new.object)
  return(new.object)
}

#' @rdname ModuleSignificance
#' @method ModuleSignificance default
#' @export
#' 
ModuleSignificance.default = function(object, ...) {
  if (inherits(object, "CorrelationNetwork")) {
    ModuleSignificance.CorrelationNetwork(object, ...)
  } else {
    stop("This method is associated with class CorrelationNetwork.")
  }
}

