## **************************************************************************
##
##    (c) 2018-2023 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    **Eigenmap class definition**
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
#' Class and Methods for Spatial Eigenvector Maps
#' 
#' Create and handle spatial eigenvector maps of a set of locations a space with
#' an arbitrary number of dimensions.
#' 
#' @name eigenmap-class
#' 
#' @docType class
#' 
#' @param x an `eigenmap-class` object.
#' @param ... Further parameters to be passed to other functions or methods
#' (currently ignored).
#' 
#' @details
#' The `print` method provides the number of the number of orthonormal
#' variables (i.e. basis functions), the number of observations these functions
#' are spanning, and their associated eigenvalues.
#' 
#' The `plot` method provides a plot of the eigenvalues and offers the
#' possibility to plot the values of variables for 1- or 2-dimensional sets of
#' coordinates. \code{plot.eigenmap} opens the default graphical device driver,
#' i.e., \code{X11}, \code{windows}, or \code{quartz} and recurses through
#' variable with a left mouse click on the graphical window. A right mouse click
#' interrupts recursing on \code{X11} and \code{windows} (Mac OS X users should
#' hit \emph{Esc} on the \code{quartz} graphical device driver (Mac OS X users).
#' 
#' @format `eigenmap-class` objects contain:
#' \describe{
#'   \item{coordinates}{ A matrix of coordinates. }
#'   \item{truncate}{ The interval within which pairs of sites are considered as
#'   neighbours. }
#'   \item{D}{ A distance matrix. }
#'   \item{weighting}{ The weighting function that had been used. }
#'   \item{wpar}{ The weighting function parameter that had been used. }
#'   \item{lambda}{ A vector of the eigenvalues obtain from the computation
#'   of the eigenvector map. }
#'   \item{U}{ A matrix of the eigenvectors defining the eigenvector map. }
#' }
#' 
#' @references
#' Borcard, D. and Legendre, P. 2002. All-scale spatial analysis of ecological
#' data by means of principal coordinates of neighbour matrices. Ecol. Model.
#' 153: 51-68
#' 
#' Dray, S.; Legendre, P. and Peres-Neto, P. 2006. Spatial modelling: a
#' comprehensive framework for principal coordinate analysis of neighbor
#' matrices (PCNM). Ecol. Modelling 196: 483-493
#' 
#' Legendre, P. and Legendre, L. 2012. Numerical Ecology, 3rd English edition.
#' Elsevier Science B.V., Amsterdam, The Netherlands.
#' 
#' @author \packageAuthor{codep}
#' Maintainer: \packageMaintainer{codep}
#' 
#' @seealso
#' \code{\link{MCA}}
#' \code{\link{eigenmap}}
#' 
#' @importFrom graphics plot
#' 
NULL
#' 
#' @describeIn eigenmap-class
#' 
#' Print method for eigenmap-class objects
#' 
#' @export
print.eigenmap <- function(x, ...) {
  cat(
    paste(
      "\nMoran's eigenvector map containing",
      length(x$lambda),
      "basis functions.\n"
    )
  )
  cat(
    paste(
      "Functions span",
      nrow(x$U),
      "observations.\n\n"
    )
  )
  cat(paste("Eigenvalues:\n"))
  print.default(x$lambda)
  cat("\n")
  return(invisible(NULL))
}
#' 
#' @describeIn eigenmap-class
#' 
#' Plot method for eigenmap-class objects
#' 
#' @export
plot.eigenmap <- function(x, ...) {
  if (ncol(x$coordinates) > 2)
    warning(
      paste(
        ncol(x$coordinates),
        "dimensions were provided but only the first 2 were used",
        "for plotting."
      )
    )
  cat(
    "Left-click on the graphical display to see further variables or",
    "right-click (Mac: esc) to terminate plotting.\n"
  )
  for (i in 1:length(x$lambda)) {
    layout(matrix(c(1,2,2),1,3))
    plot(
      y=x$lambda,
      x=1:length(x$lambda),
      ylab=expression(lambda),
      xlab="Order",
      type="b"
    )
    title("Eigenvalues diagram")
    points(
      y=x$lambda[i],
      x=i,
      pch=21L,
      bg="black"
    )
    if (ncol(x$coordinates) == 1L)
      plot(
        y=x$U[,i],
        x=x$coordinates[,1L],
        ylab="value",
        xlab="Location",
        type="l"
      )
    if (ncol(x$coordinates) > 1L) {
      plot(
        y=x$coordinates[,2L],
        x=x$coordinates[,1L],
        asp=1,
        ylab="Location (y)",
        xlab="Location (x)",
        type="p",
        pch=3L
      )
      gcol <- grey((sign(x$U[,i])+1)/2)
      gsize <- 3*abs(x$U[,i])/max(abs(x$U[,i]))
      points(
        y=x$coordinates[,2L],
        x=x$coordinates[,1L],
        pch=21L,
        bg=gcol,
        cex=gsize
      )
    }
    title(
      paste(
        "Variable display:",
        colnames(x$U)[i]
      )
    )
    ttt <- locator(1L)
    if (is.null(ttt)) break
  }
  return(invisible(NULL))
}
#' 
