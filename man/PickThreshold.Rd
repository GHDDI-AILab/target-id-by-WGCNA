% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/network.R
\name{PickThreshold}
\alias{PickThreshold}
\alias{PickThreshold.ExpAssayFrame}
\alias{PickThreshold.default}
\title{Pick a threshold.}
\usage{
PickThreshold(object, ...)

\method{PickThreshold}{ExpAssayFrame}(assay, powerVector = c(1:10, seq(12, 20, 2)))

\method{PickThreshold}{default}(object, ...)
}
\arguments{
\item{object}{An object.}

\item{...}{further arguments to be passed to or from other methods.}

\item{assay}{An ExpAssayFrame object.}

\item{powerVector}{A numeric vector of candidate soft thresholding powers for which
the scale free topology fit indices are to be calculated.}
}
\value{
A new ExpAssayFrame object.
}
\description{
Pick a threshold.

Pick a soft-thresholding power for correlation network construction.
}
