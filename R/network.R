#' @include util.R
#' @include generics.R
#' 
NULL

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

#' @rdname GetNetwork
#' @method GetNetwork ExpAssayFrame
#' @export
#' 
GetNetwork.ExpAssayFrame = function(
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
  return(net)
}

#' @rdname GetNetwork
#' @method GetNetwork default
#' @export
#'
GetNetwork.default = function(object, ...) {
  if (inherits(object, "ExpAssayFrame")) {
    GetNetwork.ExpAssayFrame(object, ...)
  } else {
    stop("This method is associated with class ExpAssayFrame.")
  }
}

