% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/preprocessing.R
\name{Reshape}
\alias{Reshape}
\alias{Reshape.ProteinGroups}
\alias{Reshape.default}
\title{Reshape a table.}
\usage{
Reshape(object, ...)

\method{Reshape}{ProteinGroups}(object)

\method{Reshape}{default}(object, ...)
}
\arguments{
\item{object}{An object of class ProteinGroups.}

\item{...}{further arguments to be passed to or from other methods.}
}
\value{
An object of class ExpAssayFrame, whose rows correspond
to samples and columns to genes.
}
\description{
Extract only the expression levels and reshape the rows and
columns for further analysis.
}
\examples{
\dontrun{
Assay = ReadProteinGroups(".")
assay = Reshape(Assay)
}

}
