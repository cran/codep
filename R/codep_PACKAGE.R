## **************************************************************************
##
##    (c) 2018-2023 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    **Package codep description**
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
#' \packageTitle{codep}
#' 
#' @description \packageDescription{codep}
#' Multiscale Codependence Analysis (MCA) consists in assessing the coherence of
#' pairs of variables in space (or time) using the product of their correlation
#' coefficients with series of spatial (or temporal) eigenfunctions. That
#' product, which is positive or negative when variables display similar or
#' opposing trends, respectively, is called a codependence coefficient.
#' 
#' The eigenfunctions used in the calculation are obtained in three steps: 1) a
#' distance matrix is calculated from the locations of samples in space (or the
#' sampling organisation through time). 2) From that distance matrix, a matrix
#' of Moran spatial weights is obtained; this is the same matrix as used to
#' calculate Moran's autocorrelation index, hence the name. And 3) the spatial
#' weight matrix is eigenvalue-decomposed after centring the rows and columns of
#' the spatial weight matrix.
#' 
#' The statistical significance of codependence coefficients is tested using
#' parametric or permutational testing of a tau statistic. The `tau` statistic
#' is the product of the Student's `t` statistics obtained from comparison of
#' the two variables with a given eigenfunction. The `tau` statistic can take
#' either positive or negative values, thereby allowing one to perform
#' one-tailed or two-tailed testing. For multiple response variables, testing is
#' performed using the `phi` statistic instead of `tau`. That statistic follows
#' the distribution of the product of two Fisher-Snedocor F statistics (see
#' \link{product-distribution} for details).
#' 
#' @name codep_PACKAGE
#' 
#' @aliases codep-package
#' 
#' @details Function \code{\link{MCA}} performs Multiscale Codependence Analysis
#' (MCA). Functions \code{\link{test.cdp}} and \code{\link{permute.cdp}} handle
#' parametric or permutation testing of the codependence coefficients,
#' respectively.
#' 
#' Methods are provided to print and plot \link{cdp-class} objects
#' (\code{\link{print.cdp}} and \code{\link{plot.cdp}}, respectively) as well as
#' summary (\code{\link{summary.cdp}}), fitted values
#' (\code{\link{fitted.cdp}}), residuals (\code{\link{residuals.cdp}}), and for
#' making predictions (\code{\link{predict.cdp}}).
#' 
#' Function \code{\link{eigenmap}} calculates spatial eigenvector maps following
#' the approach outlined in Dray et al. (2006), and which are necessary to
#' calculate \code{\link{MCA}}. It returns a \link{eigenmap-class} object. The
#' package also features methods to print (\code{\link{print.eigenmap}}) and
#' plot (\code{\link{plot.eigenmap}}) these objects. Function
#' \code{\link{eigenmap.score}} can be used to make predictions for spatial
#' models built from the eigenfunctions of \code{\link{eigenmap}} using
#' distances between one or more target locations and the sampled locations for
#' which the spatial eigenvector map was built.
#' 
#' The package also features an exemplary dataset \link{salmon} containing 76
#' sampling site positions along a 1520 m river segment. It also contains
#' functions \code{\link{cthreshold}} and \code{\link{minpermute}}, which
#' compute the testwise type I error rate threshold corresponding to a given
#' familywise threshold and the minimal number of permutations needed for
#' testing Multiscale Codependence Analysis given the alpha threshold,
#' respectively.
#' 
#' The DESCRIPTION file:
#' \packageDESCRIPTION{codep}
#' \packageIndices{codep}
#' 
#' @author \packageAuthor{codep}
#' Maintainer: \packageMaintainer{codep}
#' 
#' @references
#' Dray, S.; Legendre, P. and Peres-Neto, P. 2006. Spatial modelling: a
#' comprehensive framework for principal coordinate analysis of neighbor
#' matrices (PCNM). Ecol. Modelling 196: 483-493
#' 
#' Guénard, G., Legendre, P., Boisclair, D., and Bilodeau, M. 2010. Multiscale
#' codependence analysis: an integrated approach to analyse relationships across
#' scales. Ecology 91: 2952-2964
#' 
#' Guénard, G. Legendre, P. 2018. Bringing multivariate support to multiscale
#' codependence analysis: Assessing the drivers of community structure across
#' spatial scales. Meth. Ecol. Evol. 9: 292-304
#' 
#' @seealso
#' Legendre, P. and Legendre, L. 2012. Numerical Ecology, 3rd English edition.
#' Elsevier Science B.V., Amsterdam, The Neatherlands.
#' 
#' @examples
#' data(mite)
#' emap <- eigenmap(x = mite.geo, weighting = wf.RBF, wpar = 0.1)
#' emap
#' 
#' ## Organize the environmental variables
#' mca0 <- MCA(Y = log1p(mite.species), X = mite.env, emobj = emap)
#' mca0_partest <- test.cdp(mca0, response.tests = FALSE)
#' summary(mca0_partest)
#' plot(mca0_partest, las = 2, lwd = 2)
#' plot(mca0_partest, col = rainbow(1200)[1L:1000], las = 3, lwd = 4,
#'      main = "Codependence diagram", col.signif = "white")
#' 
#' rng <- list(x = seq(min(mite.geo[,"x"]) - 0.1, max(mite.geo[,"x"]) + 0.1, 0.05),
#'             y = seq(min(mite.geo[,"y"]) - 0.1, max(mite.geo[,"y"]) + 0.1, 0.05))
#' grid <- cbind(x = rep(rng[["x"]], length(rng[["y"]])),
#'               y = rep(rng[["y"]], each = length(rng[["x"]])))
#' newdists <- matrix(NA, nrow(grid), nrow(mite.geo))
#' for(i in 1L:nrow(grid)) {
#'   newdists[i,] <- ((mite.geo[,"x"] - grid[i,"x"])^2 +
#'                    (mite.geo[,"y"] - grid[i,"y"])^2)^0.5
#' }
#' 
#' spmeans <- colMeans(mite.species)
#' pca0 <- svd(log1p(mite.species) - rep(spmeans, each = nrow(mite.species)))
#' 
#' prd0 <- predict(
#'   mca0_partest,
#'   newdata = list(target = eigenmap.score(emap, newdists))
#' )
#' Uprd0 <- (prd0 - rep(spmeans, each = nrow(prd0))) %*% pca0$v %*%
#'   diag(pca0$d^-1)
#' 
#' ## Printing the response variable
#' prmat <- Uprd0[,1L]
#' dim(prmat) <- c(length(rng$x), length(rng$y))
#' zlim <- c(min(min(prmat), min(pca0$u[,1L])), max(max(prmat), max(pca0$u[,1L])))
#' image(z = prmat, x = rng$x, y = rng$y, asp = 1, zlim = zlim,
#'       col = rainbow(1200L)[1L:1000], ylab = "y", xlab = "x")
#' points(
#'   x = mite.geo[,"x"], y = mite.geo[,"y"], pch = 21,
#'   bg = rainbow(1200L)[round(1+(999*(pca0$u[,1L] - zlim[1L])/
#'                     (zlim[2L] - zlim[1L])),0)]
#' )
#' 
NULL
#' 
