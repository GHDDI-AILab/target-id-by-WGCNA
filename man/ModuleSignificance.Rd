% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/traits.R
\name{ModuleSignificance}
\alias{ModuleSignificance}
\alias{ModuleSignificance.CorrelationNetwork}
\alias{ModuleSignificance.default}
\title{Compute module-trait correlation and significance.}
\usage{
ModuleSignificance(object, ...)

\method{ModuleSignificance}{CorrelationNetwork}(object, samples, traits, prefix)

\method{ModuleSignificance}{default}(object, ...)
}
\arguments{
\item{object}{A CorrelationNetwork object.}

\item{...}{further arguments to be passed to or from other methods.}

\item{samples}{A character vector specifying the rows for analysis.}

\item{traits}{A character vector specifying the columns of traits for analysis.}

\item{prefix}{(A length-1 character) A prefix representing the disease.}
}
\value{
A new CorrelationNetwork object.
}
\description{
Compute module-trait correlation and significance.

Compute module-trait correlation and significance.
}
