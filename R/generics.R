#' Tidy the values.
#' 
#' @param object An object.
#' @param ... further arguments to be passed to or from other methods.
#' @export
#' 
Tidy = function(object, ...) {
  UseMethod("Tidy")
}

#' Return subsets of assays.
#' 
#' @param object An object.
#' @param ... further arguments to be passed to or from other methods.
#' @export
#' 
Subset = function(object, ...) {
  UseMethod("Subset")
}

#' Add the phenotype information to assays.
#' 
#' @param object An object.
#' @param ... further arguments to be passed to or from other methods.
#' @export
#' 
AddPhenotype = function(object, ...) {
  UseMethod("AddPhenotype")
}

#' Perform quality control for assays.
#' 
#' @param object An object.
#' @param ... further arguments to be passed to or from other methods.
#' @export
#' 
QC = function(object, ...) {
  UseMethod("QC")
}

#' Reshape a table.
#' 
#' @param object An object.
#' @param ... further arguments to be passed to or from other methods.
#' @export
#' 
Reshape = function(object, ...) {
  UseMethod("Reshape")
}

#' Perform log transformation.
#' 
#' @param object An object.
#' @param ... further arguments to be passed to or from other methods.
#' @export
#' 
LogTransform = function(object, ...) {
  UseMethod("LogTransform")
}

#' Perform normalization.
#' 
#' @param object An object.
#' @param ... further arguments to be passed to or from other methods.
#' @export
#' 
Normalize = function(object, ...) {
  UseMethod("Normalize")
}

#' Plot a distribution.
#' 
#' @param object An object.
#' @param ... further arguments to be passed to or from other methods.
#' @export
#' 
Histogram = function(object, ...) {
  UseMethod("Histogram")
}

#' Perform sample clustering.
#' 
#' @param object An object.
#' @param ... further arguments to be passed to or from other methods.
#' @export
#' 
SampleTree = function(object, ...) {
  UseMethod("SampleTree")
}

#' Pick a threshold.
#' 
#' @param object An object.
#' @param ... further arguments to be passed to or from other methods.
#' @export
#' 
PickThreshold = function(object, ...) {
  UseMethod("PickThreshold")
}

#' Construct a correlation network.
#' 
#' @param object An object.
#' @param ... further arguments to be passed to or from other methods.
#' @export
#' 
AddNetwork = function(object, ...) {
  UseMethod("AddNetwork")
}

#' Plot a gene tree.
#' 
#' @aliases GeneTree
#' @param object An object.
#' @param ... further arguments to be passed to or from other methods.
#' @export
#' 
ModulePlot = function(object, ...) {
  UseMethod("ModulePlot")
}

#' @export
#' 
GeneTree = function(...) {
  ModulePlot(...)
}

#' Compute the degree of each node.
#' 
#' @param object An object.
#' @param ... further arguments to be passed to or from other methods.
#' @export
#' 
AddConnectivity = function(object, ...) {
  UseMethod("AddConnectivity")
}

#' Return the degree of each node.
#' 
#' @param object An object.
#' @param ... further arguments to be passed to or from other methods.
#' @export
#' 
GetConnectivity = function(object, ...) {
  UseMethod("GetConnectivity")
}

#' Find hub genes.
#' 
#' @param object An object.
#' @param ... further arguments to be passed to or from other methods.
#' @export
#' 
GetHubGenes = function(object, ...) {
  UseMethod("GetHubGenes")
}

#' Compute module-trait correlation and significance.
#' 
#' @param object An object.
#' @param ... further arguments to be passed to or from other methods.
#' @export
#' 
ModuleSignificance = function(object, ...) {
  UseMethod("ModuleSignificance")
}

#' Combine module-trait correlation and significance.
#' 
#' @param object An object.
#' @param ... further arguments to be passed to or from other methods.
#' @export
#' 
BindModuleSignificance = function(object, ...) {
  UseMethod("BindModuleSignificance")
}

#' Plot a module-trait heatmap.
#' 
#' @param object An object.
#' @param ... further arguments to be passed to or from other methods.
#' @export
#' 
ModuleTraitHeatmap = function(object, ...) {
  UseMethod("ModuleTraitHeatmap")
}

#' Compute gene-trait correlation and significance.
#' 
#' @param object An object.
#' @param ... further arguments to be passed to or from other methods.
#' @export
#' 
GeneSignificance = function(object, ...) {
  UseMethod("GeneSignificance")
}

#' Compute module membership of genes.
#' 
#' @param object An object.
#' @param ... further arguments to be passed to or from other methods.
#' @export
#' 
ModuleMembership = function(object, ...) {
  UseMethod("ModuleMembership")
}

