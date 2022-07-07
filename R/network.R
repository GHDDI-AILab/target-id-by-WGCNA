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
  ATTR_POW = "powerEstimate"
  new.assay = data.table::copy(assay)
  attr(new.assay, ATTR_POW) = list()
  WGCNA::enableWGCNAThreads(nThreads = 4)
  for (i in 1:length(assay)) {
    fit = WGCNA::pickSoftThreshold(assay[[i]], powerVector = powerVector)
    attr(new.assay, ATTR_POW)[[i]] = fit$powerEstimate
  }
  names(attr(new.assay, ATTR_POW)) = names(assay)
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
    if (length(power) == 1) {
      power = rep(power, length(assay))
    } else if (length(power) != length(assay)) {
      stop("The given power was invalid!")
    }
  } else {
    stop("The given power was invalid!")
  }
  new.assay = data.table::copy(assay)
  ATTR_NET = "network"
  attr(new.assay, ATTR_NET) = list()
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
    merge = WGCNA::mergeCloseModules(assay[[i]], colors = dynamicColors, 
      cutHeight = MEDissThres)
    ## Result
    attr(new.assay, ATTR_NET)[[i]] = list(
      power = power[[i]], 
      MEDissThres = MEDissThres, 
      minModuleSize = minModuleSize, 
      adjacency = adjMat, 
      dissTOM = dissTOM, 
      geneTree = geneTree, 
      moduleEigengenes = merge$newMEs, 
      moduleColors = stats::setNames(merge$colors, genes), 
      moduleLabels = stats::setNames(match(merge$colors, colorOrder)-1, genes), 
      unmergedColors = stats::setNames(dynamicColors, genes)
      )
  }
  names(attr(new.assay, ATTR_NET)) = names(new.assay)
  class(new.assay) = c("CorrelationNetwork", class(new.assay))
  return(new.assay)
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

#' Calculate the connectivity in a network.
#' 
#' Calculate both the whole network and 
#' the intramodular connectivity for all genes.
#' 
#' @param object A CorrelationNetwork object.
#' @param geneinfo (A length-1 logical) Add gene annotation or not.
#' @return A new CorrelationNetwork object.
#' 
#' @rdname AddConnectivity
#' @method AddConnectivity CorrelationNetwork
#' @export
#' 
AddConnectivity.CorrelationNetwork = function(
  object, 
  geneinfo = FALSE
) {
  ATTR_NET = "network"
  ATTR_CON = "connectivity"
  GENE = "gene"
  #MODULE = "module"
  if (length(object) != length(attr(object, ATTR_NET))) {
    stop("Invalid CorrelationNetwork object in the input!")
  }
  new.object = data.table::copy(object)
  attr(new.object, ATTR_CON) = list()
  for (i in 1:length(new.object)) {
    network = attr(new.object, ATTR_NET)[[i]]
    module_labels = network$moduleLabels %>% 
      data.table::as.data.table(., keep.rownames = TRUE) %>% 
      data.table::setnames(., c(GENE, "module"))
  
    ## WGCNA::intramodularConnectivity() may not return the gene names 
    ## when the input includes gene names with a dash '-', like HLA-A.
    connectivity = 
      WGCNA::intramodularConnectivity(network$adjacency, network$moduleLabels) %>% 
      data.table::as.data.table(., keep.rownames = TRUE) %>% 
      data.table::setnames(., "rn", GENE) %>% 
      .[, (GENE) := module_labels[[GENE]]] %>% 
      .[module_labels, on = GENE] %>% 
      .[order(module, -kWithin)]
    
    attr(new.object, ATTR_CON)[[i]] = 
      if (geneinfo) get_geneinfo(connectivity) else connectivity
  }
  names(attr(new.object, ATTR_CON)) = names(new.object)
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

#' Return the connectivity in a network.
#' 
#' Return both the whole network and 
#' the intramodular connectivity for all genes.
#' 
#' @param object A CorrelationNetwork object.
#' @param geneinfo (A length-1 logical) Add gene annotation or not.
#' @param index A length-1 numeric or character vector specifying the frame. (default: 1)
#' @return A new CorrelationNetwork object.
#' 
#' @rdname GetConnectivity
#' @method GetConnectivity CorrelationNetwork
#' @export
#' 
GetConnectivity.CorrelationNetwork = function(
  object, 
  geneinfo = TRUE, 
  index = 1
) {
  ATTR_CON = "connectivity"
  if (length(object) != length(attr(object, ATTR_CON))) {
    object = AddConnectivity(object)
  }
  connectivity = attr(object, ATTR_CON)[[
    assert_length_1(index)
    ]]
  if (geneinfo) get_geneinfo(connectivity) else connectivity
}

#' @rdname GetConnectivity
#' @method GetConnectivity default
#' @export
#' 
GetConnectivity.default = function(object, ...) {
  if (inherits(object, "CorrelationNetwork")) {
    GetConnectivity.CorrelationNetwork(object, ...)
  } else {
    stop("This method is associated with class CorrelationNetwork.")
  }
}

#' Find hub genes in each module of a correlation network.
#' 
#' @param object A CorrelationNetwork object.
#' @param intramodular.ratio.threshold Choose the genes with the connectivity higher 
#'   than 0.9 * top_connectivity_within_the_module, 
#'   when intramodular.ratio.threshold = 0.9.
#' @param ratio.threshold Choose the genes with the connectivity higher 
#'   than 0.9 * top_connectivity_across_all_modules, 
#'   when ratio.threshold = 0.9.
#' @param top.n (A length-1 integer) Choose n top genes in each modules.
#' @param index A length-1 numeric or character vector specifying the frame. (default: 1)
#' @param geneinfo (A length-1 logical) Add gene annotation or not.
#' @return A data.table with hub genes.
#' 
#' @rdname GetHubGenes
#' @method GetHubGenes CorrelationNetwork
#' @export
#' 
GetHubGenes.CorrelationNetwork = function(
  object, 
  intramodular.ratio.threshold = 0.9, 
  ratio.threshold = 0.0, 
  top.n, 
  index = 1, 
  geneinfo = TRUE
) {
  ATTR_CON = "connectivity"
  if (length(object) != length(attr(object, ATTR_CON))) {
    object = AddConnectivity.CorrelationNetwork(object)
  }
  connectivity = attr(object, ATTR_CON)[[
    assert_length_1(index)[[1]]
    ]]
  if (missing(top.n) || length(top.n) < 1) {
    hub_genes = connectivity[, max_kWithin := max(kWithin)
      ][, top_kWithin := max(kWithin), by = "module"
      ][module  > 0 & 
        kWithin > top_kWithin*intramodular.ratio.threshold & 
        kWithin > max_kWithin*ratio.threshold
      ][, max_kWithin := NULL
      ][, top_kWithin := NULL
      ][order(module, -kWithin)]
  } else {
    top.n = assert_length_1(top.n)
    hub_genes = connectivity[order(module, -kWithin)
      ][, head(.SD, top.n), by = "module"]
  }
  if (geneinfo) get_geneinfo(hub_genes) else hub_genes
}

#' @rdname GetHubGenes
#' @method GetHubGenes default
#' @export
#' 
GetHubGenes.default = function(object, ...) {
  if (inherits(object, "CorrelationNetwork")) {
    GetHubGenes.CorrelationNetwork(object, ...)
  } else {
    stop("This method is associated with class CorrelationNetwork.")
  }
}

