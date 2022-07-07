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
  if (missing(prefix)) {
    prefix = ""
  } else {
    prefix = paste0(assert_length_1(prefix), ".")
  }
  if (missing(traits) || length(traits) < 1) {
    traits = colnames(attr(object, ATTR_PHENO))
  } else if (! all(traits %in% colnames(attr(object, ATTR_PHENO)))) {
    stop("Invalid traits in the input!")
  }
  new.object = data.table::copy(object)
  attr(new.object, ATTR_MS) = list()
  for (i in 1:length(new.object)) {
    if (missing(samples)) {
      MEs = attr(object, ATTR_NET)[[i]][["moduleEigengenes"]]
      datTraits = attr(object, ATTR_PHENO)[rownames(MEs), traits, drop = FALSE]
    } else if (length(samples) > 1) {
      MEs = attr(new.object, ATTR_NET)[[i]][["moduleEigengenes"]][samples, ]
      datTraits = attr(object, ATTR_PHENO)[samples, traits, drop = FALSE]
    } else {
      stop("The number of samples should be more than one!")
    }
    module_trait_cor = WGCNA::cor(MEs, datTraits, use = "pairwise")
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

#' Combine the module-trait correlation and significance results.
#' 
#' Combine the module-trait correlation and significance results 
#' from two CorrelationNetwork objects.
#' 
#' @param this A CorrelationNetwork object.
#' @param that Another CorrelationNetwork object.
#' @return A new CorrelationNetwork object.
#' @export
#' 
BindModuleSignificance.CorrelationNetwork = function(
  this, 
  that 
) {
  ATTR_MS = "module-trait"
  if (! inherits(that, "CorrelationNetwork")) {
    stop("Need two CorrelationNetwork objects as input.")
  }
  if (length(this) != length(that) || 
      ! setequal(names(attr(this, ATTR_MS)), names(attr(that, ATTR_MS)))) {
    stop("The module-trait of the two input objects should correspond!")
  }
  if (length(this) != length(attr(this, ATTR_MS))) {
    stop("Invalid CorrelationNetwork object: 'this'!")
  }
  if (length(that) != length(attr(that, ATTR_MS))) {
    stop("Invalid CorrelationNetwork object: 'that'!")
  }
  if (! is.null(names(attr(this, ATTR_MS))) && ! is.null(names(attr(that, ATTR_MS)))) {
    indices = names(attr(this, ATTR_MS))
  } else {
    indices = 1:length(this)
  }
  merge_module_trait = function(df1, df2) {
    d.f = merge(df1, df2, by = "row.names", all = TRUE)
    rownames(d.f) = d.f[["Row.names"]]
    d.f[["Row.names"]] = NULL
    return(d.f)
  }
  new.object = data.table::copy(this)
  for (i in indices) {
    attr(new.object, ATTR_MS)[[i]] = list(
      cor = merge_module_trait(attr(this, ATTR_MS)[[i]]$cor, attr(that, ATTR_MS)[[i]]$cor), 
      pval = merge_module_trait(attr(this, ATTR_MS)[[i]]$pval, attr(that, ATTR_MS)[[i]]$pval)
      )
  }
  names(attr(new.object, ATTR_MS)) = names(attr(this, ATTR_MS))
  return(new.object)
}

#' @rdname BindModuleSignificance
#' @method BindModuleSignificance default
#' @export
#' 
BindModuleSignificance.default = function(object, ...) {
  if (inherits(object, "CorrelationNetwork")) {
    BindModuleSignificance.CorrelationNetwork(object, ...)
  } else {
    stop("This method is associated with class CorrelationNetwork.")
  }
}

#' Get hub genes in the modules associated with traits of interest.
#' 
#' @param object A CorrelationNetwork object.
#' @param traits A character vector specifying the columns of traits for analysis.
#' @param prefix (A length-1 character) A prefix representing the disease.
#' @param index A length-1 numeric or character vector specifying the frame. (default: 1)
#' @param ... further arguments to be passed to `GetHubGenes()`.
#' @return A data.table with hub genes.
#' @export
#' 
GetRelatedHubGenes.CorrelationNetwork = function(
  object, 
  traits, 
  prefix, 
  index = 1, 
  ...
) {
  ATTR_MS = "module-trait"
  ATTR_NET = "network"
  if (length(object) != length(attr(object, ATTR_MS)) || 
      length(object) != length(attr(object, ATTR_NET))) {
    stop("Invalid CorrelationNetwork object was given!")
  }
  if (missing(prefix)) {
    prefix = ""
  } else {
    prefix = paste0(assert_length_1(prefix), ".")
  }
  pval = attr(object, ATTR_MS)[[
    assert_length_1(index)[[1]]
    ]]$pval
  if (missing(traits) || length(traits) < 1) {
    traits = colnames(pval)
  } else {
    traits = paste0("p.", prefix, traits)
  }
  nModules = nrow(pval)
  get_modules = function(col) rownames(pval)[pval[[col]] < 0.05/nModules]
  moduleColors = lapply(traits, get_modules) %>% unlist %>% sub("^ME", "", .)
  colorOrder = c("grey", WGCNA::standardColors(50))
  moduleLabels = match(moduleColors, colorOrder) - 1
  GetHubGenes(object, index = index, ...)[module %in% moduleLabels, ]
}

#' @rdname GetRelatedHubGenes
#' @method GetRelatedHubGenes default
#' @export
#' 
GetRelatedHubGenes.default = function(object, ...) {
  if (inherits(object, "CorrelationNetwork")) {
    GetRelatedHubGenes.CorrelationNetwork(object, ...)
  } else {
    stop("This method is associated with class CorrelationNetwork.")
  }
}

