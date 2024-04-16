## **************************************************************************
##
##    (c) 2018-2023 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    **Geodesic Distances Calculation**
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
#' Calculation of Geodesic Distances
#' 
#' Function \code{geodesics} carries out the calculation of pairwise geodesic
#' distances within a set of coordinates or between two sets thereof, using one
#' of two calculation approaches.
#' 
#' @param x A set of geographic coordinates in the form of a two-column
#' \code{\link{matrix}} or \code{\link{data.frame}}.
#' @param y An other two-column \code{\link{matrix}} or \code{\link{data.frame}}
#' containing an optional second set of coordinates.
#' @param method The calculation method used to obtain the distances (default:
#' haversine method; see details).
#' @param radius Radius of the planetary body (when assuming a sphere; default:
#' 6371000 m).
#' @param sma Length of the semi-major axis of the planetary body (when assuming
#' a revolution ellipsoid; default: 6378137 m).
#' @param flat Flattening of the ellipsoid (default: 1/298.257223563).
#' @param maxiter Maximum number of iterations, whenever iterative calculation
#' is involved (default: 1024).
#' @param tol Tolerance used when iterative calculation is involved (default:
#' \code{.Machine$double.eps^0.75}; a machine dependent value).
#' 
#' @return A `\link{dist}-class` object or, whenever argument \code{y} is
#' provided, a \code{\link{matrix}} with as many rows as the number of rows in
#' argument \code{x} and as many columns as the number of rows in argument
#' \code{y}.
#' 
#' @details When only one set of coordinates is given to the function (i.e.,
#' when argument \code{y} is omitted), the function returns the pairwise
#' distances in the form of a `\link{dist}-class` object representing a
#' lower-triangle matrix. When the second coordinate set is given, the function
#' calculates the distances between each coordinate of argument \code{x} and
#' each coordinate of argument \code{y}.
#' 
#' Two calculation methods are implemented. The first is the haversine formula,
#' which assume the planetary body to be a sphere. The radius of that sphere is
#' given to the function as its argument \code{radius}, with the default value
#' being the mean radius of planet earth. Of the two methods implemented, the
#' haversine formula is fastest but its precision depend on how well the
#' planetary body match the sphericity assumption. The second method implemented
#' is Vincenty's inverse formula, which assumes the the planetary body is a
#' revolution ellipsoid, which is expected for rotating semi-fluid such as
#' planet earth. Argument \code{sma}, the length of the semi-major axis,
#' corresponds to the radius of the circle obtained when the revolution
#' ellipsoid at the equator, whereas argument \code{flat} correspond to the
#' compression of the sphere, along the diameter joining the poles, to form the
#' ellipsoid of revolution. Their default values corresponds to parameters for
#' planet Earth according to WGS84. These values, along with arguments
#' \code{maxiter} and \code{tol}, are ignored when using the haversine formula,
#' whereas the value of argument \code{radius} is ignored when using Vincenty's
#' inverse formula.
#' 
#' Vincenty's inverse formula is more precise on planet Earth (on the order of
#' 0.5mm) than the haversine formula, but it involves more computation time and
#' may sometimes fail to converge. This is more likely for pairs of locations
#' that are nearly antipodal or both (numerically) very close to the equator.
#' The results returned by the function when using Vincenty's inverse formula
#' are given a \code{niter} attribute that gives the number of iterations that
#' were necessary to achieve convergence. Numbers greater than argument
#' \code{maxiter} are indicative of failed convergence; a warning is issued in
#' such a circumstance.
#' 
#' Geodesic distance matrices are non metric.
#' 
#' @author \packageAuthor{codep}
#' Maintainer: \packageMaintainer{codep}
#' 
#' @seealso The \code{\link{dist}-class} and associated methods.
#' 
#' @references
#' Vincenty, T. 1975. Direct and Inverse Solutions of Geodesics on the Ellipsoid 
#' with application of nested equations. Survey Review XXIII (176): 88-93
#' doi:10.1179/sre.1975.23.176.88
#' 
#' Inman, J. 1835. Navigation and Nautical Astronomy: For the Use of British
#' Seamen (3 ed.). London, UK: W. Woodward, C. & J. Rivington
#' 
#' @examples
#' ##
#' ### First example: locations spread throughout the world
#' ##
#' coords <- cbind(c(43,22,9,12,-40,72,-86,-22),
#'                 c(-135,22,0,1,-45,12,27,-139))
#' res_hav <- geodesics(coords)  ## Default: the haversine formula
#' res_hav
#' res_vif <- geodesics(coords, method = "Vincenty")
#' res_vif
#' attr(res_vif,"niter") ## The numbers of iterations
#' res_vif-res_hav       ## Absolute difference
#' 200*(res_vif-res_hav)/(res_vif+res_hav) ## Large relative difference
#' ##
#' ### Second example: locations nearer from one another
#' ##
#' coords <- cbind(c(45.01,44.82,45.23,44.74),
#'                 c(72.03,72.34,71.89,72.45))
#' res_hav <- geodesics(coords)
#' res_vif <- geodesics(coords, method = "Vincenty")
#' res_vif-res_hav       ## Absolute difference
#' 200*(res_vif-res_hav)/(res_vif+res_hav) ## Relative difference are smaller
#' ##
#' 
#' @useDynLib codep, .registration = TRUE
#' 
#' @export
geodesics <- function(x, y, method = c("haversine", "Vincenty"),
                      radius = 6.371e6, sma = 6378137.0, flat = 1/298.257223563,
                      maxiter = 1024L, tol = .Machine$double.eps^0.75) {
  if (NCOL(x)!=2L)
    stop("'x' must be Lat Lon geodesic coordinates!")
  if (!is.matrix(x))
    x <- as.matrix(x)
  storage.mode(x) <- "double"
  method <- match.arg(method)
  if (!missing(y)) {
    if(NCOL(y)!=2L)
      stop("'y' must be Lat Lon geodesic coordinates!")
    if(!is.matrix(y))
      y <- as.matrix(y)
    storage.mode(y) <- "double"
    N <- c(NROW(x),NROW(y))
    out <- matrix(0, N[2L], N[1L])
    rownames(out) <- dimnames(y)[[1L]]
    colnames(out) <- dimnames(x)[[1L]]
  } else {
    N <- c(NROW(x),NROW(x))
    out <- double(N[1L]*(N[1L]-1L)/2)
  }
  if (method == "haversine")
    res <- .C(
      "dist_geo_hvs",
      x,
      if(missing(y)) x else y,
      N,
      missing(y),
      out,
      radius
    )
  if (method == "Vincenty") {
    niter <- integer(length(out))
    res <- .C(
      "dist_geo_vif",
      x,
      if(missing(y)) x else y,
      N,
      missing(y),
      out,
      niter,
      sma,
      flat,
      maxiter,
      tol
    )
  }
  out[] <- res[[5L]]
  attr(out, "method") <- sprintf("geodesic:%s", method)
  attr(out, "call") <- match.call()
  if (missing(y)) {
    attr(out, "Size") <- N[1L]
    attr(out, "Labels") <- dimnames(x)[[1L]]
    attr(out, "Diag") <- FALSE
    attr(out, "Upper") <- FALSE
    if (method == "Vincenty")
      attr(out, "niter") <- res[[6L]]
    attr(out, "class") <- "dist"
  }
  if (method == "Vincenty" && any(res[[6L]] > maxiter))
    warning(
      "Convergence failed for ",
      sum(res[[6L]]>maxiter),
      " distance(s)!"
    )
  return(out)
}
#' 
