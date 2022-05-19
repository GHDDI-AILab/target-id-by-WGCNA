#' Target Identification Using the WGCNA Method
#' 
#' Target identification is an essential first step in drug
#'   discovery. This package implements convenient functions for performing
#'   target identification tasks on gene expression data using the WGCNA method.
#' 
#' Functions:
#' - `ReadExperimentalDesign()`
#' - `ReadProteinGroups()`
#' 
#' Classes:
#' - `ExperimentList-class`
#'   - `ExperimentInfo-class`
#'   - `ExpAssayTable-class`
#'     - `ProteinGroups-class`
#'   - `ExpAssayFrame-class`
#' 
#' S3 methods:
#' - `Subset()`
#' - `QC()`
#' - `Reshape()`
#' - `LogNorm()`
#' 
#' @importFrom utils str capture.output
#' @importFrom stats dist as.dist
#' @importFrom data.table data.table as.data.table fread fwrite setnames ":=" copy
#' @importFrom magrittr "%>%" "%T>%" "%<>%"
#' @importFrom WGCNA adjacency TOMsimilarity labels2colors moduleEigengenes mergeCloseModules plotDendroAndColors standardColors pickSoftThreshold intramodularConnectivity cor corPvalueStudent labeledHeatmap verboseScatterplot
#' @importFrom fastcluster hclust
#' @importFrom dynamicTreeCut cutreeDynamic
#' 
#' @docType package
#' @keywords internal
"_PACKAGE"

# Make sure data.table knows we are using it
.datatable.aware = TRUE
