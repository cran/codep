\encoding{utf8}
\name{mca-class}
\alias{mca-class}
\alias{print.mca}
\alias{summary.mca}
\alias{plot.mca}
\alias{fitted.mca}
\alias{residuals.mca}
\alias{predict.mca}
\docType{class}
\title{ Class and methods for Multiscale Codependence Analysis }
\description{
  Class and methods to handle Multiscale Codependence Analysis (MCA)
}
\usage{
\method{print}{mca}(x, ...)
\method{summary}{mca}(object, ...)
\method{plot}{mca}(x, ...)
\method{fitted}{mca}(object, which=NA, components=FALSE, ...)
\method{residuals}{mca}(object, which=NA, ...)
\method{predict}{mca}(object, which=NA, newdata=NA, components=FALSE, ...)
}
\arguments{
  \item{x}{ A \code{\link{mca-class}} object. }
  \item{object}{ A \code{\link{mca-class}} object. }
  \item{which}{ A numeric vector of indices or character vector variable
    names to test or force-use. Mandatory if \code{mcaobj} is
    untested. }
  \item{components}{ A boolean specifying whether the components of
    fitted or predicted values associated with single eigenfunctions in
    the map should be returned. }
  \item{newdata}{ A numeric vector containing new values of the
    explonatory variable. }
  \item{...}{ Further parameters to be passed to other functions or
    methods (currently ignored). }
}
\details{
  The \code{fitted}, \code{residuals}, and \code{predict} methods return
  a single-column matrix of fitted, residuals, or predicted values,
  respectively. The \code{fitted} and \code{predict} methods return a
  list a list when the parameter \code{component} is \code{TRUE}. The
  list contains the \code{fitted} or \code{predicted} values as a first
  item and a matrix \code{components} as a second. This matrix has one
  column for each statistically significant codependence coefficient.
}
\value{
  \code{\link{mca-class}} objects contain:
  \item{data}{ A copy of the response (\code{y}) and explanatory
    (\code{x}) variables that were given to \code{\link{MCA}}. }
  \item{emobj}{ The \code{\link{eigenmap-class}} object that was given to
    \code{\link{MCA}}. }
  \item{Upyxcb}{ A 4 columns matrix containing the vectors of
    cross-products of structuring variable (\code{U}) with variable
    \code{y} and \code{x}, the codependence coefficients \eqn{\mathbf{c}}{c} and
    the coregression coefficients \eqn{\mathbf{b}}{b}. }
  \item{test}{ Results of statistical testing as performed by
    \code{\link{test.mca}} or \code{\link{permute.mca}}. \code{NULL} if
    no testing was performed, such as when only \code{\link{MCA}} had
    been called. The results of statistical testing is a list containing
    the following members: }
  \item{permute}{ The number of randomized permutations used by
    \code{permute.mca} for permutation testing. \code{0} or \code{FALSE}
    for parametric testing obtained using \code{\link{test.mca}}. }
  \item{significant}{ The indices of codependence coefficient describing
    statistically significant codependence between \code{y} and
    \code{x}, in decreasing order of magnitude. }
  \item{test.table}{ The testing table (a 4 columns matrix) with
    \eqn{\tau}{tau} statistics, degrees-of-freedom, and testwise and
    familywise probabilities of type I (\eqn{\alpha}{alpha}) error. It
    contains one line for each statistically significant coefficient (if
    any). }
  \item{details}{ Details about permutation testing not shown in
    \code{test.table}. \code{NULL} for parametric testing. }

}
\references{
  Guénard, G., Legendre, P., Boisclair, D., and Bilodeau, M. 2010.
  Multiscale codependence analysis: an integrated approach to analyse
  relationships across scales. Ecology 91: 2952-2964
}
\author{ Guillaume Guénard, Département des sciences biologiques,
  Université de Montréal, Montréal, Québec, Canada. }
