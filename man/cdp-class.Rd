\encoding{utf8}
\name{cdp-class}
\alias{cdp}
\alias{cdp-class}
\alias{print.cdp}
\alias{summary.cdp}
\alias{plot.cdp}
\alias{fitted.cdp}
\alias{residuals.cdp}
\alias{predict.cdp}
\docType{class}
\title{ Class and methods for Multiscale Codependence Analysis involving
multiple descriptors }
\description{
  Class and methods to handle Multiscale Codependence Analysis (mMCA)
}
\usage{
\method{print}{cdp}(x, ...)
\method{plot}{cdp}(x, col, col.signif=2, main="", ...)
\method{summary}{cdp}(object, ...)
\method{fitted}{cdp}(object, selection, components=FALSE, ...)
\method{residuals}{cdp}(object, selection, ...)
\method{predict}{cdp}(object, selection, newdata, components=FALSE, ...)
}
\arguments{
  \item{x, object}{ A \code{\link{cdp-class}} object. }
  \item{col}{ A vector of color values to be used for plotting the
    multivariate codependence coefficients. }
  \item{col.signif}{ Color of the frame used to mark the statistically
    significant codependence coefficients . }
  \item{main}{ Text for the main title of the plot. }
  \item{selection}{ A numeric vector of indices or character vector variable
    names to test or force-use. Mandatory if \code{object} is
    untested. }
  \item{components}{ A boolean specifying whether the components of
    fitted or predicted values associated with single eigenfunctions in
    the map should be returned. }
  \item{newdata}{ A list with elements $X, $meanY, and $target that
    contain the information needed to make predictions (see details). }
  \item{...}{ Further parameters to be passed to other functions or
    methods. }
}
\details{
  The \code{fitted}, \code{residuals}, and \code{predict} methods return
  a matrix of fitted, residuals, or predicted values, respectively. The
  \code{fitted} and \code{predict} methods return a list a list when the
  parameter \code{component} is \code{TRUE}. The list contains the
  \code{fitted} or \code{predicted} values as a first element and an
  array \code{components} as a second. That 3-dimensional array has one
  matrix for each statistically significant codependence coefficient.

  For making predictions, parameter \code{newdata} may contain three
  elements: \code{$X}, a matrix of new values of the explanatory
  variables, \code{$meanY}, a vector of the predicted mean values of the
  responses, and \code{$target}, a matrix of target scores for arbitraty
  locations within the study area. When no \code{$X} is supplied, the
  descriptor given to \code{\link{MCA}} is recycled, while when no
  \code{$meanY} is supplied, the mean values of the response variables
  given to \code{\link{MCA}} are used. Finally, when element
  \code{$target} is omitted from \code{newdata}, predictions are made at
  the sites were observations were done. When none of the above is
  provided, or if \code{newdata} is omitted when calling the prediction
  method, the behaviour of the \code{predict} method is identical to
  that of the \code{fitted} method.

  From version 0.7-1, \code{\link{cdp-class}} replaces the former class
  \code{mca} used by \code{\link{codep-package}} because the standard
  package \code{MASS} also had S3 methods for a class named \code{mca}
  that were overwritten by those of \code{\link{codep-package}}.
}
\value{
  \code{\link{cdp-class}} objects contain:
  \item{data}{ A list with two elements: the first being a copy of the
    response (\code{Y}) and the second being a copy of the explanatory
    variables (\code{X}). This is the variables that were given to
    \code{\link{MCA}}. }
  \item{emobj}{ The \code{\link{eigenmap-class}} object that was given to
    \code{\link{MCA}}. }
  \item{UpYXcb}{ A list with five elements: the first (\code{UpY}) is a
    matrix of the cross-products of structuring variable (\code{U}) and
    the response variable \code{Y}, the second (\code{UpX}) is a matrix
    of the cross-product of the structuring variable and the explanatory
    variables (\code{X}), the third (\code{C}) is a 3-dimensional array
    of the codependence coefficients, the fourth (\code{B}) is a
    3-dimensional array of the coregression coefficients, and the fifth
    (\code{CM}) is a matrix of the multivariate codependence
    coefficients. }
  \item{test}{ Results of statistical testing as performed by
    \code{\link{test.cdp}} or \code{\link{permute.cdp}}. \code{NULL}
    if no testing was performed, such as when only \code{\link{MCA}}
    had been called. The results of statistical testing is a list
    containing the following members: }
  \item{$permute}{ The number of randomized permutations used by
    \code{permute.cdp} for permutation testing. \code{0} or \code{FALSE}
    for parametric testing obtained using \code{\link{test.cdp}}. }
  \item{$significant}{ The indices of codependence coefficient describing
    statistically significant codependence between \code{Y} and
    \code{X}, in decreasing order of magnitude. }
  \item{$global}{ The testing table (a 5-column matrix) with
    \eqn{\phi}{phi} statistics, degrees-of-freedom, and testwise and
    familywise probabilities of type I (\eqn{\alpha}{alpha}) error. It
    contains one line for each statistically significant global
    coefficient (if any) in addition to test results for the first,
    non-significant coefficient, on which the testing procedure
    stopped. }
  \item{$response}{ Tests of every single response variable (a
    3-dimensional array), had such tests been requested while calling
    the testing function, \code{NULL} otherwise. }
  \item{$permutations}{ Details about permutation testing not shown in
    \code{test$global} or \code{test$response}. \code{NULL} for parametric
    testing. }
}
\references{
  Guénard, G., Legendre, P., Boisclair, D., and Bilodeau, M. 2010.
  Multiscale codependence analysis: an integrated approach to analyse
  relationships across scales. Ecology 91: 2952-2964
  
  Guénard, G. Legendre, P. 2018. Bringing multivariate support to
  multiscale codependence analysis: Assessing the drivers of community
  structure across spatial scales. Meth. Ecol. Evol. 9: 292-304
}
\author{ Guillaume Guénard, Département des sciences biologiques,
  Université de Montréal, Montréal, Québec, Canada. }
