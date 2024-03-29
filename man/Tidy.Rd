% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/preprocessing.R
\name{Tidy}
\alias{Tidy}
\alias{Tidy.ExpAssayTable}
\alias{Tidy.default}
\title{Tidy the values.}
\usage{
Tidy(object, ...)

\method{Tidy}{ExpAssayTable}(object)

\method{Tidy}{default}(object, ...)
}
\arguments{
\item{object}{An object of class ExpAssayTable.}

\item{...}{further arguments to be passed to or from other methods.}
}
\value{
A new object of class ExpAssayTable.
}
\description{
Tidy the values.

Tidy an ExpAssayTable object.
}
\examples{
\dontrun{
new.Assay = Tidy(old.Assay)
}

}
