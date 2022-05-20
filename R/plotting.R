#' @rdname DistributionPlot
#' @method DistributionPlot ExpAssayFrame
#' @export
#' 
DistributionPlot.ExpAssayFrame = function(
  assay, 
  samples, 
  file
) {
  if (missing(file)) {
    file = format(Sys.Date(), "Plots_%Y%m%d/LevelDistribution.pdf")
  }
  if (!missing(samples) && length(samples)) {
    assay = Subset(assay, samples)
  }
  assay = assert_length_1(assay)
  d.f = data.frame(unlist(assay[[1]], use.names = FALSE))
  names(d.f) = "Level"
  p = ggplot2::ggplot(d.f, ggplot2::aes(x = Level)) +   
    ggplot2::geom_histogram(ggplot2::aes(y = ..count..), binwidth = 0.05, 
                            alpha = .3, position = "identity")
  dir.create(dirname(file), showWarnings = FALSE)
  ggplot2::ggsave(file, p)
}

#' @rdname DistributionPlot
#' @method DistributionPlot default
#' @export
#' 
DistributionPlot.default = function(object, ...) {
  if (inherits(object, "ExpAssayFrame")) {
    DistributionPlot.ExpAssayFrame(object, ...)
  } else {
    stop("This method is associated with class ExpAssayFrame.")
  }
}

#' Make a plot of sample clustering using an expression profile.
#' 
#' @param assay An ExpAssayFrame object.
#' @param file A length-1 character naming the file to draw the plot.
#' @rdname SampleTree
#' @method SampleTree ExpAssayFrame
#' @export
#' 
SampleTree.ExpAssayFrame = function(
  assay, 
  file, 
  title = "Sample Clustering", 
  hline = 50, 
  plot.width = 20, 
  plot.height = 12
) {
  tree = assay %>% 
    assert_length_1() %>% 
    .[[1]] %>% 
    stats::dist() %>% 
    fastcluster::hclust(., method = "average")
  if (missing(file)) {
    file = format(Sys.Date(), "Plots_%Y%m%d/SampleClustering.pdf")
  }
  dir.create(dirname(file), showWarnings = FALSE)
  pdf(file = file, width = plot.width, height = plot.height)
  par(cex = 0.6)
  par(mar = c(5, 10, 10, 0))
  plot(tree, main = title, sub = "", xlab = "", 
       cex.lab = 2, cex.axis = 2, cex.main = 4)
  abline(h = hline, col = "red")
  dev.off()
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

#' Plot the gene tree.
#' 
#' @rdname ModulePlot
#' @method ModulePlot ExpAssayFrame
#' @export
#' 
ModulePlot.ExpAssayFrame = function(
  assay, 
  file.prefix, 
  plot.width = 12.5, 
  plot.height = 10
) {
  if (missing(file.prefix)) {
    file.prefix = ""
  } else if (is.character(file.prefix)) {
    file.prefix = paste0(assert_length_1(file.prefix), ".")
  } else {
    stop("The given file.prefix was invalid!")
  }
  net = attr(assay, "Network") %>% assert_length_1() %>% .[[1]]
  dir = format(Sys.Date(), "Plots_%Y%m%d")
  filename = sprintf("%s/%sGeneTree.power%s.pdf", dir, file.prefix, net$power)
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  pdf(filename, width = plot.width, height = plot.height)
  WGCNA::plotDendroAndColors(
    dendro = net$geneTree, colors = cbind(net$unmergedColors, net$moduleColors), 
    groupLabels = c("Dynamic Tree Cut", "Module colors"),
    dendroLabels = FALSE, hang = 0.03,
    addGuide = TRUE, guideHang = 0.05, 
    cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5,
    main = "Gene dendrogram and module colors"
    )
  dev.off()
  invisible(filename)
}

#' @rdname ModulePlot
#' @method ModulePlot default
#' @export
#' 
ModulePlot.default = function(object, ...) {
  if (inherits(object, "ExpAssayFrame")) {
    ModulePlot.ExpAssayFrame(object, ...)
  } else {
    stop("This method is associated with class ExpAssayFrame.")
  }
}

