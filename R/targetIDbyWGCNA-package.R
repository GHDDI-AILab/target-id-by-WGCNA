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
#' @importFrom graphics abline par
#' @importFrom grDevices dev.off pdf
#' @importFrom utils capture.output head str
#' @importFrom stats as.dist dist setNames
#' @importFrom data.table ":=" .SD as.data.table copy data.table fread fwrite setnames
#' @importFrom magrittr "%>%" "%T>%"
#' @importFrom ggplot2 aes geom_histogram ggplot ggsave
#' @importFrom WGCNA TOMsimilarity adjacency cor corPvalueStudent intramodularConnectivity labeledHeatmap labels2colors moduleEigengenes mergeCloseModules plotDendroAndColors pickSoftThreshold standardColors verboseScatterplot
#' @importFrom fastcluster hclust
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom AnnotationDbi select
#' 
#' @docType package
#' @keywords internal
"_PACKAGE"

# Make sure data.table knows we are using it
.datatable.aware = TRUE
