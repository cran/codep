## **************************************************************************
##
##    (c) 2018-2023 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    **Weighting Functions for Spatial Eigenvector Maps**
##
##    This file is part of codep
##
##    codep is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    codep is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with codep. If not, see <https://www.gnu.org/licenses/>.
##
##    R source code file
##
## **************************************************************************
##
#' Weighting Functions for Spatial Eigenvector Map
#' 
#' A set of common distance weighting functions to calculate spatial eignevector
#' maps using function \code{\link{eigenmap}}.
#' 
#' @name weighting-functions
#' 
#' @param d A triangular (`\link{dist}-class`) or rectangular geographic
#' distance matrix produced by \code{\link{dist}}, \code{\link{Euclid}}, or
#' \code{\link{geodesics}}.
#' @param boundaries Where applicable, a two-element numeric vector containing
#' the lower and upper threshold values used to obtain the connectivity matrix.
#' (see details).
#' @param wpar Where applicable, a parameter controlling the shape of the
#' spatial weighting function.
#' 
#' @returns A `\link{dist}-class` object when argument \code{d} is a
#' `\link{dist}-class` object or a rectangular matrix when argument \code{d} is
#' a rectangular matrix, either one with the weights as its values.
#' 
#' @details These functions are meant primarily to be called within functions
#' \code{\link{eigenmap}} and \code{\link{eigenmap.score}}. In
#' \code{\link{eigenmap}}, argument \code{d} is a lower-triangular
#' `\link{dist}-class` object and the resulting lower-triangular weight matrix
#' is used in calculating the spatial eigenvector map. In
#' \code{\link{eigenmap.score}}, \code{d} is a rectangular matrix of the
#' distances between a set of arbitrary locations (rows) and reference locations
#' (columns; the locations for which the the spatial eigenvector map has been
#' built and the resulting rectangular weight matrix is used to calculate
#' spatial eigenfunction values. These values allow one to use the spatial
#' information of a data set for making predictions at arbitrary values.
#' 
#' `Wf.sqrd` (default value) consists in taking w[i,j] = -0.5*d[i,j] and does
#' not involve any truncation.
#' 
#' `Wf.RBF` consists in taking \code{w[i,j] = exp(-wpar*d[i,j]^2)} and
#' does not involve any truncation, where \code{wpar} is a non-zero real
#' positive value (default: 1).
#' 
#' `Wf.binary` the spatial weighting matrix is simply the connectivity
#' matrix.
#' 
#' `Wf.PCNM` is \code{a[i,j] = 1 - (d[i,j] / (wpar*boundaries[2]))^2}, where
#' \code{wpar} is a non-zero real positive value (default: 4).
#'  
#' `Wf.Drayf1` is \code{a[i,j] = 1 - (d[i,j] / dmax)} where \code{d_max} is the
#' distance between the two most distant locations in the set.
#' 
#' `Wf.Drayf2` is \code{a[i,j] = 1 - (d[i,j] / dmax)^wpar}, where \code{wpar} is
#' a non-zero real positive value (default: 1).
#' 
#' `Wf.Drayf3` is \code{a[i,j] = 1 / d[i,j]^wpar}, where \code{wpar} is a
#' non-zero real positive value (default: 1).
#' 
#' Functions \code{Wf.Drayf1}, \code{Wf.Drayf2}, and \code{Wf.Drayf3} were
#' proposed by Dray et al. (2006) and function \code{PCNM} was proposed by
#' Legendre and Legendre (2012).
#' 
#' The \code{Wf.sqrd} weighting approach is equivalent to submitting the
#' elementwise square-root of the distance matrix to a principal coordinate
#' analysis. It was proposed by Diniz-Filho et al. (2013) and is equivalent, for
#' evenly spaced transect or surfaces (square or rectangle), to using the basis
#' functions of type II discrete cosine basis transforms; a fact that has gone
#' unnoticed by Diniz-Filho et al. (2013).
#' 
#' The radial basis function (RBF) is a widespread kernel method involving sets
#' of real-valued functions whose values depend on the distance between any
#' given input coordinate and a set of fixed points (a single fixed point for
#' each function). It is implemented using function \code{Wf.RBF} using all the
#' sampling points as the fixed points.
#' 
#' When calculating the connectivity matrix, pairs of location whose distance to
#' one another are between the boundary values (argument \code{bounraries}) are
#' considered as neighbours (\code{b[i,j]=1}) whereas values located below the
#' minimum and above the maximum are considered as equivalent or distant,
#' respectively (\code{b[i,j]=0} in both cases).
#' 
#' User may implement custom weighting functions. These functions must at the
#' very least have an argument \code{d}, and can be given arguments
#' \code{boundaries} and \code{wpar}. Argument \code{wpar} may be a vector with
#' any number of elements. They should be added to the R-code file
#' (weighting-functions.R). User-provided weighting functions with an argument
#' \code{wpar} must come with a valid default value for that parameter since
#' \code{\link{eigenmap}} may internally call it without a formal value.
#' 
#' @author \packageAuthor{codep}
#' Maintainer: \packageMaintainer{codep}
#' 
#' @references
#' Borcard, D. and Legendre, P. 2002. All-scale spatial analysis of ecological
#' data by means of principal coordinates of neighbour matrices. Ecol. Model.
#' 153: 51-68
#' 
#' Diniz-Filho, J. A. F.; Diniz, J. V. B. P. L.; Rangel, T. F.; Soares, T. F.;
#' de Campos Telles, M. P.; Garcia Collevatti, R. and Bini, L. M. 2013. A new
#' eigenfunction spatial analysis describing population genetic structure.
#' Genetica 141:479-489.
#' 
#' Dray, S.; Legendre, P. and Peres-Neto, P. 2006. Spatial modelling: a
#' comprehensive framework for principal coordinate analysis of neighbor
#' matrices (PCNM). Ecol. Modelling 196: 483-493
#' 
#' Legendre, P. and Legendre, L. 2012. Numerical Ecology, 3rd English edition.
#' Elsevier Science B.V., Amsterdam, The Netherlands.
#' 
#' @examples
#' locations <- c(1,2,4,7,10,14,17,21)
#' D <- dist(locations)
#' wf.sqrd(D)
#' wf.RBF(D, wpar = 0.1)
#' wf.binary(D, c(0,5))
#' wf.PCNM(D, c(0,5))
#' wf.Drayf1(D, c(0,5))
#' wf.Drayf2(D, c(0,5), 0.5)
#' wf.Drayf3(D, c(0,5), 0.5)
#' 
#' emap <- eigenmap(D, locations, wf.Drayf2, c(0,5), 0.5)
#' emap
#' 
#' emap <- eigenmap(D, locations, wf.Drayf3, c(0,5), 0.25)
#' emap
#' 
#' emap <- eigenmap(D, locations, wf.RBF, wpar = 0.1)
#' emap
#' 
NULL
#' 
#' @describeIn weighting-functions
#' 
#' Principal coordinates of the square-root distance matrix (Diniz-Filho et al.
#' 2013).
#' 
#' @export
wf.sqrd <- function(d) {
  w <- -0.5 * d
  attr(w, "method") <- sprintf("sqrd(%s)", attr(w,"method"))
  attr(w, "call") <- c(attr(d, "call"), match.call())
  return(w)
}
#' 
#' @describeIn weighting-functions
#' 
#' Radial basis functions with the observations as the kernels.
#' 
#' @export
wf.RBF <- function(d, wpar = 1) {
  w <- exp(-wpar * d**2)
  w[d==0] <- 0
  attr(w, "method") <- sprintf("RBF(%s)", attr(w, "method"))
  attr(w, "call") <- c(attr(d, "call"), match.call())
  return(w)
}
#' 
#' @describeIn weighting-functions
#' 
#' Borcard & Legendre's (2002) principal coordinates of the neighbour matrix
#' approach.
#' 
#' @export
wf.PCNM <- function(d, boundaries, wpar = 4) {
  w <- d
  w[!(d > boundaries[1L] & d <= boundaries[2L])] <- wpar * boundaries[2L]
  w <- -0.5 * w**2
  attr(w, "method") <- sprintf("PCNM(%s)", attr(w, "method"))
  attr(w, "call") <- c(attr(d, "call"), match.call())
  return(w)
}
#' 
#' @describeIn weighting-functions
#' 
#' Dray et al. (2006) Moran's eigenvector maps (distance-based binary
#' connections without continuous weighting of the neighbours).
#' 
#' @export
wf.binary <- function(d, boundaries) {
  b <- d
  b[] <- 0
  b[d > boundaries[1L] & d <= boundaries[2L]] <- 1
  attr(b, "method") <- sprintf("binary(%s)", attr(d, "method"))
  attr(b, "call") <- c(attr(d, "call"), match.call())
  return(b)
}
#' 
#' @describeIn weighting-functions
#' 
#' Dray et al. (2006) Moran's eigenvector maps (distance-based binary
#' connections with continuous weighting of the neighbours: f1).
#' 
#' @export
wf.Drayf1 <- function(d, boundaries) {
  b <- d
  b[] <- 0
  b[d > boundaries[1L] & d <= boundaries[2L]] <- 1
  a <- 1 - d / max(d)
  w <- b * a
  w[d==0] <- 0
  attr(w, "method") <- sprintf("Drayf1(%s)", attr(w, "method"))
  attr(w, "call") <- c(attr(d, "call"), match.call())
  return(w)
}
#' 
#' @describeIn weighting-functions
#' 
#' Dray et al. (2006) Moran's eigenvector maps (distance-based binary
#' connections with continuous weighting of the neighbours: f2).
#' 
#' @export
wf.Drayf2 <- function(d, boundaries, wpar = 1) {
  b <- d
  b[] <- 0
  b[d > boundaries[1L] & d <= boundaries[2L]] <- 1
  a <- 1 - (d / max(d))^wpar
  w <- b * a
  w[d==0] <- 0
  attr(w, "method") <- sprintf("Drayf2(%s)", attr(w, "method"))
  attr(w, "call") <- c(attr(d, "call"), match.call())
  return(w)
}
#' 
#' @describeIn weighting-functions
#' 
#' Dray et al. (2006) Moran's eigenvector maps (distance-based binary
#' connections with continuous weighting of the neighbours: f3).
#' 
#' @export
wf.Drayf3 <- function(d, boundaries, wpar = 1) {
  b <- d
  b[] <- 0
  b[d > boundaries[1L] & d <= boundaries[2L]] <- 1
  a <- 1 / d^wpar
  w <- b * a
  w[d==0] <- 0
  attr(w, "method") <- sprintf("Drayf3(%s)", attr(w, "method"))
  attr(w, "call") <- c(attr(d, "call"), match.call())
  return(w)
}
#' 
