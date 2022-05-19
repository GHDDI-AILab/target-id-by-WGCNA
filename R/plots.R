#' Make a plot of sample clustering using an expression profile.
#' 
#' @param object An ExpAssayFrame object.
#' @param file A length-1 character naming the file to draw the plot.
#' @rdname SampleTree
#' @method SampleTree ExpAssayFrame
#' @export
#' 
SampleTree.ExpAssayFrame = function(
  object, 
  file
) {

}

#' @rdname SampleTree
#' @method SampleTree default
#' @export
SampleTree.default = function(object, ...) {
  if (inherits(object, "ExpAssayFrame")) {
    SampleTree.ExpAssayFrame(object, ...)
  } else {
    stop("This method is associated with class ExpAssayFrame.")
  }
}

