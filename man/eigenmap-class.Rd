\encoding{utf8}
\name{eigenmap-class}
\alias{eigenmap-class}
\alias{print.eigenmap}
\alias{plot.eigenmap}
\docType{class}
\title{Class and methods for spatial eigenvector maps}
\description{
  Methods to handle spatial eigenvector maps of a set of locations in a
  space with an arbitrary number of dimension.
}
\usage{
\method{print}{eigenmap}(x, ...)
\method{plot}{eigenmap}(x, ...)
}
\arguments{
  \item{x}{ An object of \code{\link{eigenmap-class}} }
  \item{...}{ Further parameters to be passed to other functions or
    methods (currently ignored). }
}
\details{
  The \code{print} method provides the number of the number of
  orthonormal variables (i.e. basis functions), the number of
  observations these functions are spanning, and their associated
  eigenvalues.
  
  The \code{plot} method provides a plot of the eigenvalues and offers
  the possibility to plot the values of variables for 1- or
  2-dimensional sets of coordinates. \code{plot.eigenmap} opens the default
  graphical device driver, i.e., \code{X11}, \code{windows}, or
  \code{quartz} and recurses through variable with a left mouse click on
  the graphical window. A right mouse click interrupts recursing on
  \code{X11} and \code{windows} (Mac OS X users should hit \emph{Esc} on
  the \code{quartz} graphical device driver (Mac OS X users).
}
\value{
  \code{\link{eigenmap-class}} objects contain:
  \item{coordinates}{ a matrix of coordinates, }
  \item{truncate}{ the interval within which pairs of sites are
  considered as neighbour, }
  \item{D}{ the distance matrix, }
  \item{weighting}{ the weighting function that had been used, }
  \item{wpar}{ the weighting function parameter that had been used, }
  \item{lambda}{ a vector of the eigenvalues obtain from the computation
    of the eigenvector map, and }
  \item{U}{ a matrix of the eigenvectors defining the eigenvector map. }
}
\references{
  Dray, S.; Legendre, P. and Peres-Neto, P. 2006. Spatial modelling: a
  comprehensive framework for principal coordinate analysis of neighbor
  matrices (PCNM). Ecol. Modelling 196: 483-493
}
\author{ Guillaume Guénard, Département des sciences biologiques,
  Université de Montréal, Montréal, Québec, Canada. }
\seealso{
  \code{\link{MCA}}
  \code{\link{eigenmap}}
}
