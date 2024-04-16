## **************************************************************************
##
##    (c) 2018-2023 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    **Euclidean distance between two point sets A and B**
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
#' Calculation of the Euclidean Distance
#' 
#' Function \code{Euclid} carries out the calculation of pairwise Euclidean
#' distances within a set of coordinates or between two sets thereof, with
#' optional weights.
#' 
#' @param x A set of coordinates in the form of a \code{\link{matrix}} or
#' \code{\link{data.frame}}.
#' @param y An optional second set of coordinates in the same dimensions as
#' argument \code{x}.
#' @param squared Should the squared Euclidean distances be returned (default:
#' FALSE).
#' 
#' @return A `\link{dist}-class` object or, whenever \code{y} is provided,
#' a \code{\link{matrix}} with as many rows as the number of rows in \code{x}
#' and as many columns as the number of rows in \code{y}.
#' 
#' @details When only one set of coordinates is given to the function (i.e.,
#' when argument \code{y} is omitted), the function returns the pairwise
#' distances in the form of a `\link{dist}-class` object representing a
#' lower-triangle matrix. If weights are omitted, the result is identical to
#' that produced by function \link{dist} with argument
#' \code{method = "euclidean"} (the function's default).
#' 
#' The standard `R` function used to calculate the Euclidean distance
#' (\code{\link{dist}}), only allows one to calculate pairwise distances between
#' the rows of a single matrix of Cartesian coordinates and return a
#' `\link{dist}-class` object, which is a one-dimensional array meant to be
#' interpreted as a lower-triangular matrix. Function \code{Euclid} can also be
#' provided two data matrices (arguments \code{x} and \code{y}) and output a
#' rectangular matrix of the Euclidean distances.
#' 
#' @author \packageAuthor{codep}
#' Maintainer: \packageMaintainer{codep}
#' 
#' @seealso The `\link{dist}-class` and associated methods.
#' 
#' @importFrom stats dist
#' 
#' @examples
#' ## A set of reference points:
#' x <- cbind(c(1,4,5,2,8,4), c(3,6,7,1,3,2))
#' dimnames(x) <- list(LETTERS[1:6], c("x", "y"))
#' 
#' ## The pairwise Euclidean distances among the reference points: 
#' d1 <- Euclid(x)
#' d1
#' 
#' ## That result is the same as that obtained from function dist:
#' d2 <- dist(x, method = "euclidean")
#' all(d1 == d2)
#' 
#' ## A second set of points:
#' y <- cbind(c(3,5,7), c(3,6,8))
#' dimnames(y) <- list(LETTERS[7:9], c("x", "y"))
#' 
#' ## The distances between the points in y (rows) and x (columns):
#' Euclid(x, y)
#' 
#' @useDynLib codep, .registration = TRUE
#' 
#' @export
Euclid <- function(x, y, squared = FALSE) {
  m <- NCOL(x)
  if (!is.matrix(x))
    x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (!missing(y)) {
    if(NCOL(y)!=m)
      stop("'y' must have the same number of coordinates as 'x'!")
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
  res <- .C(
    "dist_Euclid",
    x,
    if(missing(y)) x else y,
    N,
    missing(y),
    m,
    out,
    squared
  )
  out[] <- res[[6L]]
  attr(out, "method") <- "euclidean"
  attr(out, "call") <- match.call()
  if (missing(y)) {
    attr(out, "Size") <- N[1L]
    attr(out, "Labels") <- dimnames(x)[[1L]]
    attr(out, "Diag") <- FALSE
    attr(out, "Upper") <- FALSE
    attr(out, "class") <- "dist"
  }
  return(out)
}
#'
