#' @include util.R
#' @include generics.R
#' 
NULL

#' Pick a soft-thresholding power for correlation network construction.
#' 
#' @param assay An ExpAssayFrame object.
#' @param powerVector A numeric vector of candidate soft thresholding powers for which 
#'   the scale free topology fit indices are to be calculated.
#' @return A new ExpAssayFrame object.
#' 
#' @rdname PickThreshold
#' @method PickThreshold ExpAssayFrame
#' @export
#' 
PickThreshold.ExpAssayFrame = function(
  assay, 
  powerVector = c(1:10, seq(12, 20, 2))
) {
  new.assay = data.table::copy(assay)
  attr(new.assay, "powerEstimate") = list()
  for (name in names(assay)) {
    fit = WGCNA::pickSoftThreshold(assay[[name]], powerVector = powerVector)
    attr(new.assay, "powerEstimate")[[name]] = fit$powerEstimate
  }
  return(new.assay)
}

#' @rdname PickThreshold
#' @method PickThreshold default
#' @export
#' 
PickThreshold.default = function(object, ...) {
  if (inherits(object, "ExpAssayFrame")) {
    PickThreshold.ExpAssayFrame(object, ...)
  } else {
    stop("This method is associated with class ExpAssayFrame.")
  }
}

#' Construct a correlation network by the gene expression profiles.
#' 
#' @param assay An ExpAssayFrame object.
#' @param power (A length-1 numeric) Soft thresholding power 
#'   for the adjacency function in WGCNA analysis.
#' @return A CorrelationNetwork object.
#' 
#' @rdname AddNetwork
#' @method AddNetwork ExpAssayFrame
#' @export
#' 
AddNetwork.ExpAssayFrame = function(
  assay, 
  power
) {
  ## Usually, the threshold for cut height depends on 
  ## the circumstances (by viewing the clustering plot).
  MEDissThres = 0.15
  ## We like large modules, so we set the minimum module size relatively high:
  minModuleSize = 30
  ## Check the given power
  if (missing(power)) {
    if (length(attr(assay, "powerEstimate")) != length(assay)) {
      assay = PickThreshold(assay)
    }
    power = attr(assay, "powerEstimate")
  } else if (is.numeric(power)) {
    power = rep(power, length(assay))
  } else {
    stop("The given power was invalid!")
  }
  net = data.table::copy(assay)
  attr(net, "Network") = list()
  colorOrder = c("grey", WGCNA::standardColors(50))
  for (i in 1:length(assay)) {

  genes = colnames(assay[[i]])
  adjMat = WGCNA::adjacency(assay[[i]], power = power[[i]])
  dissTOM = 1 - WGCNA::TOMsimilarity(adjMat)
  geneTree = fastcluster::hclust(stats::as.dist(dissTOM), method = "average")
  dynamicModules = dynamicTreeCut::cutreeDynamic(
    dendro = geneTree, distM = dissTOM, 
    deepSplit = 2, pamRespectsDendro = FALSE, 
    minClusterSize = minModuleSize
    )
  dynamicColors = WGCNA::labels2colors(dynamicModules)
  MEList = WGCNA::moduleEigengenes(assay[[i]], colors = dynamicColors)
  merge = WGCNA::mergeCloseModules(assay[[i]], colors = dynamicColors, 
    cutHeight = MEDissThres)
  ## Result
  attr(net, "Network")[[i]] = list(
    power = power[[i]], 
    MEDissThres = MEDissThres, 
    minModuleSize = minModuleSize, 
    adjacency = adjMat, 
    dissTOM = dissTOM, 
    geneTree = geneTree, 
    MEs = merge$newMEs, 
    moduleColors = stats::setNames(merge$colors, genes), 
    moduleLabels = stats::setNames(match(merge$colors, colorOrder)-1, genes), 
    unmergedColors = stats::setNames(dynamicColors, genes)
    )

  }
  names(attr(net, "Network")) = names(net)
  class(net) = c("CorrelationNetwork", class(net))
  return(net)
}

#' @rdname AddNetwork
#' @method AddNetwork default
#' @export
#'
AddNetwork.default = function(object, ...) {
  if (inherits(object, "ExpAssayFrame")) {
    AddNetwork.ExpAssayFrame(object, ...)
  } else {
    stop("This method is associated with class ExpAssayFrame.")
  }
}

#' Calculate the connectivity values for genes.
#' 
#' @param object A CorrelationNetwork object.
#' @return A new CorrelationNetwork object.
#' 
#' @rdname AddConnectivity
#' @method AddConnectivity CorrelationNetwork
#' @export
#' 
AddConnectivity.CorrelationNetwork = function(
  object
) {
  GENE = "gene"
  #MODULE = "module"
  new.object = data.table::copy(object)
  attr(new.object, "connectivity") = list()
  for (i in 1:length(new.object)) {
    network = attr(new.object, "Network")[[i]]
    module_labels = network$moduleLabels %>% 
      data.table::as.data.table(., keep.rownames = TRUE) %>% 
      data.table::setnames(., c(GENE, "module"))
    attr(new.object, "connectivity")[[i]] = 
      WGCNA::intramodularConnectivity(network$adjacency, network$moduleLabels) %>% 
      data.table::as.data.table(., keep.rownames = TRUE) %>% 
      data.table::setnames(., "rn", GENE) %>% 
      .[module_labels, on = c(GENE)] %>% 
      .[order(module, -kWithin)] %>% 
      get_geneinfo()
  }
  names(attr(new.object, "connectivity")) = names(new.object)
  return(new.object)
}

#' @rdname AddConnectivity
#' @method AddConnectivity default
#' @export
#' 
AddConnectivity.default = function(object, ...) {
  if (inherits(object, "CorrelationNetwork")) {
    AddConnectivity.CorrelationNetwork(object, ...)
  } else {
    stop("This method is associated with class CorrelationNetwork.")
  }
}

