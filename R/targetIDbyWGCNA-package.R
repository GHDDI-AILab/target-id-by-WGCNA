#' @title Target Identification Using the WGCNA Method
#' @name targetIDbyWGCNA
#' @description Target identification is an essential first step in drug
#'   discovery. This package implements convenient functions for performing
#'   target identification tasks on gene expression data using the WGCNA method.
#' 
#' @importFrom utils str
#' @importFrom data.table fread fwrite setnames
#' @importFrom magrittr "%>%" "%T>%" "%<>%"
#' 
#' @docType package
#' @aliases targetIDbyWGCNA targetIDbyWGCNA-package
#' @keywords internal
"_PACKAGE"

# Make sure data.table knows we know we're using it
.datatable.aware = TRUE
