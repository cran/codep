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
#' Spatial Eigenvector Maps
#' 
#' Function to calculate spatial eigenvector maps of a set of locations
#' in a space with an arbitrary number of dimension.
#' 
#' @name eigenmap
#' 
#' @param x A set of coordinates defined in one (numeric vector) or many (a
#' coordinate x dimension matrix) dimensions or, alternatively, a distance
#' matrix provided by \code{\link{dist}}. 
#' @param alt.coord Coordinates to be used when a distance matrix is
#' provided as x. Used for plotting purposes.
#' @param weighting The function to obtain the edge weighting matrix (see
#' details).
#' @param boundaries When required by argument \code{weighting}, a two-element
#' numeric vector containing the lower and upper threshold values used to obtain
#' the connectivity matrix (see \link{weighting-functions}).
#' @param wpar Shape parameter for argument \code{weignting} (optional).
#' @param tol The smallest absolute eigenvalue for a spatial eigenfunctions to
#' be considered as a suitable predictor. Default:
#' \code{.Machine$double.eps^0.5} (a machine-dependent value).
#' @param emap An \link{eigenmap-class} object.
#' @param target A (generally rectangular) distance matrix between a set of
#' target locations for which spatially-explicit predictions are being made
#' (rows), and the reference locations given to function \code{eigenmap}
#' (columns). See example 2.
#' 
#' @details When function \code{eigenmap} is given coordinates as its argument
#' \code{x}, they are treated as Cartesian coordinates and the distances between
#' them are assumed to be Euclidean. Otherwise (e.g., when geodesic distances
#' are used), distances have to be provided as the argument \code{x} and
#' plotting coordinates have to be supplied as argument \code{alt.coord}.
#' 
#' The weighting function (see \link{weighting-functions}) must have the
#' distances as its first argument, optionally an argument named
#' \code{boundaries} giving the boundaries within which locations are regarded
#' as neighbours and/or an argument \code{wpar} containing any other weighting
#' function parameters.
#' 
#' Default values for argument \code{boundaries} are 0 for the minimum value and
#' \code{NA} for the maximum. For weighting functions with an argument
#' \code{bounraries}, The upper value \code{NA} indicates the function to take
#' the minimum value that allow every locations to form a single cluster
#' following single linkage clustering as a maximum value (obtained internally
#' from a call to \code{\link{hclust}}.
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
#' @examples ### Example 1: A linear transect.
#' data(salmon)
#' 
#' ## A warning is issued when no boundaries are provided for a function that
#' ## requires them.
#' ## Example:
#' map <- eigenmap(x = salmon[,"Position"], weighting = wf.binary)
#' map
#' ## plot(map)
#' 
#' ## In the following examples, boundaries are provided; they are needed by the
#' ## functions.
#' map <- eigenmap(x = salmon[,"Position"], weighting = wf.binary,
#'                 boundaries = c(0,20))
#' map
#' ## plot(map)
#' 
#' map <- eigenmap(x = salmon[,"Position"], weighting = wf.Drayf1,
#'                 boundaries = c(0,20))
#' map
#' ## plot(map)
#' 
#' map <- eigenmap(x = salmon[,"Position"], weighting = wf.Drayf2,
#'                 boundaries = c(0,20))
#' map
#' ## plot(map)
#' 
#' map <- eigenmap(x = salmon[,"Position"], weighting = wf.Drayf3,
#'                 boundaries = c(0,20), wpar = 2)
#' map
#' ## plot(map)
#' 
#' map <- eigenmap(x = salmon[,"Position"], weighting = wf.PCNM,
#'                 boundaries = c(0,20))
#' map
#' ## plot(map)
#' 
#' map <- eigenmap(x = salmon[,"Position"], weighting = wf.sqrd)
#' map
#' ## plot(map)
#' 
#' map <- eigenmap(x = salmon[,"Position"], weighting = wf.RBF, wpar = 0.001)
#' map
#' ## plot(map)
#' 
#' ### Example 2: Using predictor scores
#' 
#' smpl <- c(4,7,10,14,34,56,61,64)  # A sample to be discarded
#' map <- eigenmap(x = salmon[-smpl,"Position"], weighting = wf.sqrd)
#' scr <- eigenmap.score(
#'          map, target = as.matrix(dist(salmon[,"Position"]))[,-smpl]
#'        )
#' ## Scores of sampling points are the eigenvectors
#' scr[smpl,]
#' 
#' wh <- 5L   # You can try with other vectors.
#' plot(map$U[,wh] ~ salmon[-smpl,"Position"], ylab = expression(U[5]),
#'      xlab = "Position along transect")
#' points(y = scr[smpl,wh], x = salmon[smpl,"Position"], pch = 21L,
#'        bg = "black")
#' 
#' map <- eigenmap(x = salmon[-smpl,"Position"], weighting = wf.binary,
#'                 boundaries = c(0,20))
#' scr <- eigenmap.score(
#'          map, target = as.matrix(dist(salmon[,"Position"]))[smpl,-smpl])
#' 
#' ## Plot the 8 prediction sites along particular eigenvectors, here
#' ## eigenvector #1:
#' wh <- 1L   # One could try the other vectors.
#' plot(map$U[,wh] ~ salmon[-smpl,"Position"], ylab = expression(U[1L]),
#'      xlab = "Position along transect (m)")
#' points(y = scr[,wh], x = salmon[smpl,"Position"], pch=21L, bg = "black")
#' 
#' map <- eigenmap(x = salmon[-smpl,"Position"], weighting = wf.PCNM,
#'                 boundaries = c(0,100))
#' scr <- eigenmap.score(
#'          map, target = as.matrix(dist(salmon[,"Position"]))[smpl,-smpl]
#'        )
#' 
#' wh <- 1L   # You can try with other vectors.
#' plot(map$U[,wh] ~ salmon[-smpl,"Position"], ylab = expression(U[1]),
#'      xlab = "Position along transect (m)")
#' points(y = scr[,wh], x = salmon[smpl,"Position"], pch = 21L, bg = "black")
#' 
#' ### Example 3: A unevenly sampled surface.
#' 
#' data(mite)
#' 
#' ## Example using the principal coordinates of the square root of the
#' ## (Euclidean) distances:
#' map <- eigenmap(x = as.matrix(mite.geo), weighting = wf.sqrd)
#' map
#' ## plot(map)
#' 
#' ## Example using the radial basis functions (RBF):
#' map <- eigenmap(x = as.matrix(mite.geo), weighting = wf.RBF)
#' map
#' ## plot(map)
#' 
#' @importFrom grDevices grey
#' @importFrom graphics layout locator points title
#' @importFrom stats hclust
#' 
NULL
#' 
#' @describeIn eigenmap
#' 
#' Main function for generating an eigenmap-class object from Cartesian
#' coordinates or pairwise distances.
#' 
#' @export
eigenmap <- function(x, alt.coord = NA, weighting = wf.sqrd, boundaries, wpar,
                     tol = .Machine$double.eps^0.5) {
  ##
  if (!is.function(weighting))
    stop("Parameter 'weighting' must be a weighting function.")
  if (!is.numeric(x)) stop("Parameter 'x' must be numeric!")
  if (inherits(x, "dist")) {
    D <- as.matrix(x)
    if(all(is.na(alt.coord))) {
      alt.coord <- as.matrix(x=1:nrow(D))
      rownames(alt.coord) <- rownames(D)
    } else {
      alt.coord <- as.matrix(alt.coord)
      if(nrow(alt.coord) != nrow(D)) {
        stop(
          paste(
            "You provided",
            nrow(alt.coord),
            "alternate coordinates to reference",
            nrow(D),
            "observations!"
          )
        )
      } else
        rownames(alt.coord) <- rownames(D)
    }
  } else {
    alt.coord <- as.matrix(x)
    x <- dist(x, method="euclidean")
    D <- as.matrix(x)
    rownames(alt.coord) <- rownames(D)
  }
  wformals <- names(formals(weighting))
  if (any(wformals == "boundaries") && missing(boundaries)) {
    boundaries <- c(0,max(hclust(x, method="single")$height))
    warning(
      "No boundaries given, they were set to (",
      boundaries[1L], "; ",
      boundaries[2L],")"
    )
  }
  if (any(wformals == "boundaries")) {
    if (any(wformals == "wpar") && !missing(wpar))
      W <- weighting(D, boundaries, wpar)
    else
      W <- weighting(D, boundaries)
  } else {
    if (any(wformals == "wpar") && !missing(wpar))
      W <- weighting(D, wpar)
    else
      W <- weighting(D)
  }
  n <- nrow(W)
  colWbar <- matrix(colMeans(W), 1L, n)
  MW <- mean(W)
  term <- matrix(1,n,1L) %*% colWbar
  O <- W - term - t(term) + MW
  eigO <- eigen(O)
  variables <- eigO$vectors[,abs(eigO$values) >= tol]
  rownames(variables) <- rownames(D)
  colnames(variables) <- paste("dbMEM", 1L:ncol(variables), sep="")
  return(
    structure(
      list(
        coordinates = alt.coord,
        D = D,
        weighting = weighting,
        boundaries = if(!missing(boundaries)) boundaries else NULL,
        wpar = if(!missing(wpar)) wpar else NULL,
        W = W,
        colWbar = colWbar,
        MW = MW,
        lambda = eigO$values[abs(eigO$values) >= tol],
        U = variables
      ),
      class="eigenmap"
    )
  )
}
#' 
#' @describeIn eigenmap
#' 
#' Generate scores for arbitrary locations within the scope of an existing
#' eigenvector map.
#' 
#' @export
eigenmap.score <- function(emap, target) {
  if (!is.matrix(target))
    attr(target,"dim") <- c(1L, length(target))
  n <- ncol(emap$D)
  if(n != ncol(target))
    stop(
      "Mismatch distances to target: ",
      n,
      " but ",
      ncol(emap$D),
      " are expected."
    )
  ntgt <- nrow(target)
  #
  if (!is.null(emap$boundaries)) {
    if (!is.null(emap$wpar))
      Wnk <- emap$weighting(target, emap$boundaries, emap$wpar)
    else
      Wnk <- emap$weighting(target, emap$boundaries)
  } else {
    if (!is.null(emap$wpar))
      Wnk <- emap$weighting(target, emap$wpar)
    else
      Wnk <- emap$weighting(target)
  }
  Wnk[target==0] <- 0
  scores <- (Wnk - matrix(rowMeans(Wnk),ntgt,1L) %*% matrix(1,1L,n) -
               matrix(1,ntgt,1L) %*% emap$colWbar + emap$MW) %*%
    emap$U %*% diag(emap$lambda**(-1))
  colnames(scores) <- colnames(emap$U)
  return(scores)
}
#' 
