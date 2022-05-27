#' The ExperimentList class.
#' 
#' A list of experiment data.
#' 
#' @docType class
#' @name ExperimentList-class
#' @rdname ExperimentList-class
#' 
NULL

#' The ExperimentInfo class.
#' 
#' The experiment information of assay data.
#' 
#' @docType class
#' @name ExperimentInfo-class
#' @rdname ExperimentInfo-class
#' 
NULL

#' The ExpAssayTable class.
#' 
#' A list of data.tables of experiment data, whose 
#'   rows correspond to genes.
#' 
#' @docType class
#' @name ExpAssayTable-class
#' @rdname ExpAssayTable-class
#' 
NULL

#' The ProteinGroups class.
#' 
#' The protein-level mass spectrometry data.
#' 
#' @docType class
#' @name ProteinGroups-class
#' @rdname ProteinGroups-class
#' 
NULL

#' The ExpAssayFrame class.
#' 
#' A list of data frames of experiment data, whose 
#'   rows correspond to samples and columns to genes.
#' 
#' @docType class
#' @name ExpAssayFrame-class
#' @rdname ExpAssayFrame-class
#' 
NULL

#' The CorrelationNetwork class.
#' 
#' A list of data frames of experiment data, with 
#'   a list of corresponding correlation networks.
#' 
#' @docType class
#' @name CorrelationNetwork-class
#' @rdname CorrelationNetwork-class
#' 
NULL

#' @rdname print
#' @method print ExperimentList
#' @export
#'
print.ExperimentList = function(x) {
  cat(sprintf('An object of class %s\n\n', class(x)[1]))
  utils::str(x)
  invisible(x)
}

#' @rdname print
#' @method print ExpAssayTable
#' @export
#'
print.ExpAssayTable = function(x) {
  ## Print the class
  cat(sprintf('An object of class %s\n\n', class(x)[1]))
  ## Print the attribute 'experiments'
  samples = attr(x, 'experiments')
  cat(sprintf('%d Experiment(s): ', length(samples)))
  if (length(samples) > 3) {
    cat(sprintf('"%s", "%s", "%s", ...\n', 
                samples[1], samples[2], samples[3]))
  } else {
    cat(sprintf('"%s"\n', paste(samples, collapse = '", "')))
  }
  ## Print the dimensions of each assay
  cat(sprintf('%d Assay(s): "%s"\n', 
              length(x), paste(names(x), collapse = '", "')))
  for (i in 1:length(x)) {
    cat(sprintf('\t%d features across %d samples within assay %d.\n',
                nrow(x[[i]]), length(attr(x, 'experiments')), i))
  }
  invisible(x)
}

#' @rdname print
#' @method print ExpAssayFrame
#' @export
#'
print.ExpAssayFrame = function(x) {
  ## Print the class
  cat(sprintf('An object of class %s\n\n', class(x)[1]))
  ## Print the attribute 'experiments'
  samples = attr(x, 'experiments')
  cat(sprintf('%d Experiment(s): ', length(samples)))
  if (length(samples) > 3) {
    cat(sprintf('"%s", "%s", "%s", ...\n', 
                samples[1], samples[2], samples[3]))
  } else {
    cat(sprintf('"%s"\n', paste(samples, collapse = '", "')))
  }
  ## Print the dimensions of each assay
  cat(sprintf('%d Assay(s): "%s"\n', 
              length(x), paste(names(x), collapse = '", "')))
  for (i in 1:length(x)) {
    cat(sprintf('\t%d features across %d samples within assay %d.\n',
                ncol(x[[i]]), length(attr(x, 'experiments')), i))
  }
  ## Print other attributes
  cat("Attributes:\n")
  index = setdiff(names(attributes(x)), c("names", "experiments", "class"))
  utils::str(attributes(x)[index])
  invisible(x)
}

