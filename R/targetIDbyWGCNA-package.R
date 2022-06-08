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
#'     - `CorrelationNetwork-class`
#' 
#' S3 methods:
#' - `Tidy()`
#' - `Subset()`
#' - `QC()`
#' - `Reshape()`
#' - `LogTransform()`
#' - `Normalize()`
#' - `Histogram()`
#' - `SampleTree()`
#' - `PickThreshold()`
#' - `AddNetwork()`
#' - `ModulePlot()`
#' - `AddConnectivity()`
#' - `GetHubGenes()`
#' 
#' @importFrom graphics abline par
#' @importFrom grDevices dev.off pdf
#' @importFrom utils capture.output head str
#' @importFrom stats as.dist dist median setNames
#' @importFrom data.table ":=" .SD as.data.table copy data.table fread fwrite setnames
#' @importFrom magrittr "%>%" "%T>%"
#' @importFrom ggplot2 aes geom_histogram ggplot ggsave
#' @importFrom WGCNA TOMsimilarity adjacency intramodularConnectivity mergeCloseModules moduleEigengenes pickSoftThreshold
#' @importFrom WGCNA cor corPvalueStudent labeledHeatmap labels2colors plotDendroAndColors standardColors verboseScatterplot
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
