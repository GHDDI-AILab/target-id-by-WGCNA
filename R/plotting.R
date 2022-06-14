#' @include util.R
#' @include generics.R
#' 
NULL

#' Plot the distribution of expression levels in the assay.
#' 
#' @param object An ExpAssayFrame object.
#' @param samples A character vector specifying the samples of interest.
#' @param genes A character vector specifying the genes of interest.
#' @param index A length-1 numeric or character vector specifying the frame. (default: 1)
#' @param file (A length-1 character vector) File path to save the plot.
#' @param title (A length-1 character vector) Title for the plot.
#' @param plot.width (A length-1 numeric) Set the plot width. 
#'   If not supplied, use the size of current graphics device.
#' @param plot.height (A length-1 numeric) Set the plot height.
#'   If not supplied, use the size of current graphics device.
#' @return Path to the output file.
#' 
#' @rdname Histogram
#' @method Histogram ExpAssayFrame
#' @export
#' 
Histogram.ExpAssayFrame = function(
  object, 
  samples, 
  genes, 
  index = 1, 
  file, 
  title, 
  plot.width = NA, 
  plot.height = NA
) {
  if (missing(file)) {
    file = format(Sys.Date(), "Plots_%Y%m%d/LevelDistribution.pdf")
  } else if (is.character(file) && length(file) > 0) {
    file = file[[1]]
  } else {
    stop("The given file was invalid!")
  }
  if (missing(title)) {
    title = "Distribution of expression levels"
  } else if (is.character(title) && length(title) > 0) {
    title = title[[1]]
  } else {
    stop("The given title was invalid!")
  }
  assay = object[[
    assert_length_1(index)[[1]]
    ]]
  if (!missing(samples) && length(samples) > 0 && !any(is.na(samples))) {
    assay = assay[samples, ]
  }
  if (!missing(genes) && length(genes) > 0 && !any(is.na(genes))) {
    assay = assay[, genes]
  }
  d.f = data.frame(unlist(assay, use.names = FALSE))
  names(d.f) = "Level"
  p = ggplot2::ggplot(d.f, ggplot2::aes(x = Level)) + 
    ggplot2::geom_histogram(ggplot2::aes(y = ..count..), 
                            binwidth = 0.2, 
                            color = "black", 
                            fill = "white", 
                            #alpha = 0.3, 
                            position = "identity") + 
    ggplot2::labs(title = title) + 
    ggplot2::theme_classic(base_size = 20, base_family = "Times") + 
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(file, p, width = plot.width, height = plot.height)
  invisible(file)
}

#' @rdname Histogram
#' @method Histogram default
#' @export
#' 
Histogram.default = function(object, ...) {
  if (inherits(object, "ExpAssayFrame")) {
    Histogram.ExpAssayFrame(object, ...)
  } else {
    stop("This method is associated with class ExpAssayFrame.")
  }
}

#' Make a plot of sample clustering using an expression profile.
#' 
#' @param object An ExpAssayFrame object.
#' @param index A length-1 numeric or character vector specifying the frame. (default: 1)
#' @param file (A length-1 character vector) File path to save the plot.
#' @param title (A length-1 character vector) Title for the plot.
#' @param hline (A length-1 numeric) Set the position of a red horizontal line.
#' @param plot.width (A length-1 numeric) Set the plot width.
#' @param plot.height (A length-1 numeric) Set the plot height.
#' @return Path to the output file.
#' 
#' @rdname SampleTree
#' @method SampleTree ExpAssayFrame
#' @export
#' 
SampleTree.ExpAssayFrame = function(
  object, 
  index = 1, 
  file, 
  title = "Sample Clustering", 
  hline = 50, 
  plot.width = 20, 
  plot.height = 12
) {
  if (missing(file)) {
    file = format(Sys.Date(), "Plots_%Y%m%d/SampleClustering.pdf")
  } else if (is.character(file) && length(file) > 0) {
    file = file[[1]]
  } else {
    stop("The given file was invalid!")
  }
  d.f = object[[
    assert_length_1(index)[[1]]
    ]]
  tree = fastcluster::hclust(stats::dist(d.f), method = "average")
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  pdf(file = file, width = plot.width, height = plot.height)
  par(mar = c(3, 8, 8, 0))
  plot(tree, main = title, sub = "", xlab = "", 
       cex = 0.6, cex.axis = 2.5, cex.lab = 2.5, cex.main = 2.5)
  abline(h = hline, col = "red")
  dev.off()
  invisible(file)
}

#' @rdname SampleTree
#' @method SampleTree default
#' @export
#' 
SampleTree.default = function(object, ...) {
  if (inherits(object, "ExpAssayFrame")) {
    SampleTree.ExpAssayFrame(object, ...)
  } else {
    stop("This method is associated with class ExpAssayFrame.")
  }
}

#' Plot the gene tree of a gene co-expression network.
#' 
#' @param object A CorrelationNetwork object.
#' @param index A length-1 numeric or character vector specifying the frame. (default: 1)
#' @param file.prefix (A length-1 character) A prefix for the name of output file.
#' @param title (A length-1 character vector) Title for the plot.
#' @param plot.width (A length-1 numeric) Set the plot width.
#' @param plot.height (A length-1 numeric) Set the plot height.
#' @return Path to the output file.
#' 
#' @rdname ModulePlot
#' @method ModulePlot CorrelationNetwork
#' @export
#' 
ModulePlot.CorrelationNetwork = function(
  object, 
  index = 1, 
  file.prefix, 
  title = "Gene dendrogram and module colors", 
  plot.width = 12.5, 
  plot.height = 10
) {
  ATTR_NET = "network"
  net = attr(object, ATTR_NET)[[
    assert_length_1(index)[[1]]
    ]]
  if (missing(file.prefix)) {
    file.prefix = ""
  } else if (is.character(file.prefix)) {
    file.prefix = paste0(assert_length_1(file.prefix), ".")
  } else {
    stop("The given file.prefix was invalid!")
  }
  dir = format(Sys.Date(), "Plots_%Y%m%d")
  file = sprintf("%s/%sGeneTree.power%s.pdf", dir, file.prefix, net$power)
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  pdf(file = file, width = plot.width, height = plot.height)
  WGCNA::plotDendroAndColors(
    dendro = net$geneTree, 
    colors = cbind(net$unmergedColors, net$moduleColors), 
    groupLabels = c("Dynamic Tree Cut", "Module colors"),
    dendroLabels = FALSE, hang = 0.03,
    addGuide = TRUE, guideHang = 0.05, 
    cex.colorLabels = 1.8, cex.dendroLabels = 1, 
    cex.main = 2, cex.lab = 2, cex.axis = 2, 
    marAll = c(3, 12, 5, 3), 
    main = title[[1]]
    )
  dev.off()
  invisible(file)
}

#' @rdname ModulePlot
#' @method ModulePlot default
#' @export
#' 
ModulePlot.default = function(object, ...) {
  if (inherits(object, "CorrelationNetwork")) {
    ModulePlot.CorrelationNetwork(object, ...)
  } else {
    stop("This method is associated with class CorrelationNetwork.")
  }
}

