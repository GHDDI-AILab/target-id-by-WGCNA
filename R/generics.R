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

