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

#' Calculate the degree of each node.
#' 
#' @param object An object.
#' @param ... further arguments to be passed to or from other methods.
#' @export
#' 
AddConnectivity = function(object, ...) {
  UseMethod("AddConnectivity")
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

