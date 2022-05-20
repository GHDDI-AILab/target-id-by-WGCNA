#' Subset assays.
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

#' Perform log normalization.
#' 
#' @param object An object.
#' @param ... further arguments to be passed to or from other methods.
#' @export
#' 
LogNorm = function(object, ...) {
  UseMethod("LogNorm")
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

#' Draw the distribution.
#' 
#' @param object An object.
#' @param ... further arguments to be passed to or from other methods.
#' @export
#' 
DistributionPlot = function(object, ...) {
  UseMethod("DistributionPlot")
}

#' Pick a soft threshold.
#' 
#' @param object An object.
#' @param ... further arguments to be passed to or from other methods.
#' @export
#' 
PickThreshold = function(object, ...) {
  UseMethod("PickThreshold")
}

#' Construct a correlation network..
#' 
#' @param object An object.
#' @param ... further arguments to be passed to or from other methods.
#' @export
#' 
GetNetwork = function(object, ...) {
  UseMethod("GetNetwork")
}

#' Plot a gene tree..
#' 
#' @param object An object.
#' @param ... further arguments to be passed to or from other methods.
#' @export
#' 
ModulePlot = function(object, ...) {
  UseMethod("ModulePlot")
}

