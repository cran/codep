## **************************************************************************
##
##    (c) 2018-2023 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    **Multiscale codependence functions**
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
#' Multiple-descriptors, Multiscale Codependence Analysis
#' 
#' Class, Functions, and methods to perform Multiscale Codependence Analysis
#' (MCA)
#' 
#' @name MCA
#' 
#' @param Y A numeric matrix or vector containing the response variable(s).
#' @param X A numeric matrix or vector containing the explanatory variable(s).
#' @param emobj A \link{eigenmap-class} object.
#' @param object A \link{cdp-class} object.
#' @param alpha The type I (alpha) error threshold used by the
#' testing procedure.
#' @param max.step The maximum number of steps to perform when testing for
#' statistical significance.
#' @param response.tests A boolean specifying whether to test individual
#' response variables.
#' @param permute The number of random permutations used for testing. When
#' omitted, the number of permutations is calculated using function
#' \code{\link{minpermute}}.
#' @param nnode The number of parallel computation nodes.
#' @param seeds Seeds for computation nodes' random number generators when using
#' parallel computation during the permutation test.
#' @param verbose Whether to return user notifications.
#' @param ... Parameters to be passed to function \code{parallel::makeCluster}
#' 
#' @details Multiscale Codependence Analysis (MCA) allows to calculate
#' correlation-like (i.e.codependence) coefficients between two variables with
#' respect to structuring variables (Moran's eigenvector maps). The purpose of
#' this function is limited to parameter fitting.
#' 
#' Test procedures are handled through \code{test.cdp} (parametric testing) or
#' \code{permute.cdp} (permutation testing). Moreover, methods are provided for
#' printing (\code{print.cdp}), displaying a summary of the tests
#' (\code{summary.cdp}), plotting results (\code{plot.cdp}), calculating
#' fitted (\code{fitted.cdp}) and residuals values (\code{redisuals.cdp}), and
#' making predictions (\code{predict.cdp}).
#' 
#' It is noteworthy that the test procedure used by \code{MCA} deviates from the
#' standard R workflow since intermediate testing functions (\code{test.cdp} and
#' \code{permute.cdp}) need first to be called before any testing be performed.
#' 
#' Function \code{parPermute.cdp} allows the user to spread the number of
#' permutation on many computation nodes. It relies on package parallel.
#' Omitting argument \code{nnode} lets function \code{parallel::detectCores}
#' specify the number of node. Similarly, omitting parameter \code{seeds} lets
#' the function define the seeds as a set of values drawn from a uniform random
#' distribution between with minimum value \code{-.Machine$integer.max} and
#' maximum value \code{.Machine$integer.max}.
#' 
#' @returns A \link{cdp-class} object.
#' 
#' @references
#' Guénard, G., Legendre, P., Boisclair, D., and Bilodeau, M. 2010. Multiscale
#' codependence analysis: an integrated approach to analyse relationships across
#' scales. Ecology 91: 2952-2964
#' 
#' Guénard, G. Legendre, P. 2018. Bringing multivariate support to multiscale
#' codependence analysis: Assessing the drivers of community structure across
#' spatial scales. Meth. Ecol. Evol. 9: 292-304
#' 
#' @author \packageAuthor{codep}
#' Maintainer: \packageMaintainer{codep}
#' 
#' @examples ### Example 1: St. Marguerite River Salmon Transect
#' data(salmon)
#' 
#' ## Converting the data from data frames to to matrices:
#' Abundance <- log1p(as.matrix(salmon[,"Abundance",drop = FALSE]))
#' Environ <- as.matrix(salmon[,3L:5])
#' 
#' ## Creating a spatial eigenvector map:
#' map1 <- eigenmap(x = salmon[,"Position"], weighting = wf.binary,
#'                  boundaries = c(0,20))
#' 
#' ## Case of a single descriptor:
#' mca1 <- MCA(Y = Abundance, X = Environ[,"Substrate",drop = FALSE],
#'             emobj = map1)
#' mca1
#' mca1_partest <- test.cdp(mca1)
#' mca1_partest
#' summary(mca1_partest)
#' par(mar = c(6,4,2,4))
#' plot(mca1_partest, las = 3, lwd=2)
#' mca1_pertest <- permute.cdp(mca1)
#' \dontrun{
#' ## or:
#' mca1_pertest <- parPermute.cdp(mca1, permute = 999999)
#' }
#' mca1_pertest
#' summary(mca1_pertest)
#' plot(mca1_pertest, las = 3)
#' mca1_pertest$UpYXcb$C # Array containing the codependence coefficients
#' 
#' ## With all descriptors at once:
#' mca2 <- MCA(Y = log1p(as.matrix(salmon[,"Abundance",drop = FALSE])),
#'             X = as.matrix(salmon[,3L:5]), emobj = map1)
#' mca2
#' mca2_partest <- test.cdp(mca2)
#' mca2_partest
#' summary(mca2_partest)
#' par(mar = c(6,4,2,4))
#' plot(mca2_partest, las = 3, lwd=2)
#' mca2_pertest <- permute.cdp(mca2)
#' \dontrun{
#' ## or:
#'     mca2_pertest <- parPermute.cdp(mca2, permute = 999999)
#' }
#' mca2_pertest
#' summary(mca2_pertest)
#' plot(mca2_pertest, las = 3, lwd=2)
#' mca2_pertest$UpYXcb$C # Array containing the codependence coefficients
#' mca2_pertest$UpYXcb$C[,1L,] # now turned into a matrix.
#' 
#' ### Example 2: Doubs Fish Community Transect
#' 
#' data(Doubs)
#' 
#' ## Sites with no fish observed are excluded:
#' excl <- which(rowSums(Doubs.fish) == 0)
#' 
#' ## Creating a spatial eigenvector map:
#' map2 <- eigenmap(x = Doubs.geo[-excl,"DFS"])
#' ## The eigenvalues are in map2$lambda, the MEM eigenvectors in matrix map2$U
#' 
#' ## MCA with multivariate response data analyzed on the basis of the Hellinger
#' ## distance:
#' Y <- LGTransforms(Doubs.fish[-excl,],"Hellinger")
#' 
#' mca3 <- MCA(Y = Y, X=Doubs.env[-excl,], emobj = map2)
#' mca3_pertest <- permute.cdp(mca3)
#' \dontrun{
#' ## or:
#' mca3_pertest <- parPermute.cdp(mca3, permute = 999999)
#' }
#' 
#' mca3_pertest
#' summary(mca3_pertest)
#' par(mar = c(6,4,2,4))
#' plot(mca3_pertest, las = 2, lwd=2)
#' 
#' ## Array containing all the codependence coefficients:
#' mca3_pertest$UpYXcb$C
#' 
#' ## Display the results along the transect
#' spmeans <- colMeans(Y)
#' pca1 <- svd(Y - rep(spmeans, each=nrow(Y)))
#' par(mar = c(5,5,2,5) + 0.1)
#' plot(y = pca1$u[,1L], x = Doubs.geo[-excl,"DFS"], pch = 21L, bg = "red",
#'      ylab = "PCA1 loadings", xlab = "Distance from river source (km)")
#' 
#' ## A regular transect of sites from 0 to 450 (km) spaced by 1 km intervals
#' ## (451 sites in total). It is used for plotting spatially-explicit
#' ## predictions.
#' 
#' x <- seq(0,450,1)
#' newdists <- matrix(NA, length(x), nrow(Doubs.geo[-excl,]))
#' for(i in 1L:nrow(newdists))
#'   newdists[i,] <- abs(Doubs.geo[-excl,"DFS"] - x[i])
#' 
#' ## Calculating predictions for the regular transect under the same set of
#' ## environmental conditions from which the codependence model was built.
#' prd1 <- predict(mca3_pertest,
#'                 newdata = list(target = eigenmap.score(map2, newdists)))
#' 
#' ## Projection of the predicted species abundance on pca1:
#' Uprd1 <-
#'   (prd1 - rep(spmeans, each = nrow(prd1))) %*%
#'   pca1$v %*% diag(pca1$d^-1)
#' lines(y = Uprd1[,1L], x = x, col=2, lty = 1)
#' 
#' ## Projection of the predicted species abundance on pca2:
#' plot(y = pca1$u[,2L], x = Doubs.geo[-excl,"DFS"], pch = 21L, bg = "red",
#'      ylab = "PCA2 loadings", xlab = "Distance from river source (km)")
#' lines(y = Uprd1[,2L], x = x, col = 2, lty = 1)
#' 
#' ## Displaying only the observed and predicted abundance for Brown Trout.
#' par(new = TRUE)
#' plot(y = Y[,"TRU"], Doubs.geo[-excl,"DFS"], pch = 21L,
#'      bg = "green", ylab = "", xlab = "", new = FALSE, axes = FALSE)
#' axis(4)
#' lines(y = prd1[,"TRU"], x = x, col = 3)
#' mtext(side = 4, "sqrt(Brown trout rel. abundance)", line = 2.5)
#' 
#' ### Example 3: Borcard et al. Oribatid Mite
#' 
#' ## Testing the (2-dimensional) spatial codependence between the Oribatid Mite
#' ## community structure and environmental variables, while displaying the
#' ## total effects of the significant variables on the community structure
#' ## (i.e., its first principal component).
#' 
#' data(mite)
#' 
#' map3 <- eigenmap(x = mite.geo)
#' 
#' Y <- LGTransforms(mite.species, "Hellinger")
#' 
#' ## Organize the environmental variables
#' mca4 <- MCA(Y = Y, X = mite.env, emobj = map3)
#' mca4_partest <- test.cdp(mca4, response.tests = FALSE)
#' summary(mca4_partest)
#' plot(mca4_partest, las = 2, lwd = 2)
#' plot(mca4_partest, col = rainbow(1200)[1L:1000], las = 3, lwd = 4,
#'      main = "Codependence diagram", col.signif = "white")
#' 
#' ## Making a regular point grid for plotting the spatially-explicit
#' ## predictions:
#' rng <- list(
#'   x = seq(min(mite.geo[,"x"]) - 0.1, max(mite.geo[,"x"]) + 0.1, 0.05),
#'   y = seq(min(mite.geo[,"y"]) - 0.1, max(mite.geo[,"y"]) + 0.1, 0.05))
#' grid <- cbind(x = rep(rng[["x"]], length(rng[["y"]])),
#'               y = rep(rng[["y"]], each = length(rng[["x"]])))
#' newdists <- matrix(NA, nrow(grid), nrow(mite.geo))
#' for(i in 1L:nrow(grid)) {
#'   newdists[i,] <- ((mite.geo[,"x"] - grid[i,"x"])^2 +
#'                      (mite.geo[,"y"] - grid[i,"y"])^2)^0.5
#' }
#' 
#' spmeans <- colMeans(Y)
#' pca2 <- svd(Y - rep(spmeans, each = nrow(Y)))
#' 
#' prd2 <- predict(mca4_partest,
#'                 newdata = list(target = eigenmap.score(map3, newdists)))
#' Uprd2 <-
#'   (prd2 - rep(spmeans, each = nrow(prd2))) %*%
#'   pca2$v %*% diag(pca2$d^-1)
#' 
#' ## Printing the response variable (first principal component of the mite
#' ## community structure).
#' prmat <- Uprd2[,1L]
#' dim(prmat) <- c(length(rng$x), length(rng$y))
#' zlim <- c(min(min(prmat), min(pca2$u[,1L])), max(max(prmat),
#'                                                  max(pca2$u[,1L])))
#' image(z = prmat, x = rng$x, y = rng$y, asp = 1, zlim = zlim,
#'       col = rainbow(1200L)[1L:1000], ylab = "y", xlab = "x")
#' points(
#'   x=mite.geo[,"x"], y=mite.geo[,"y"], pch=21L,
#'   bg = rainbow(1200L)[round(1+(999*(pca2$u[,1L]-zlim[1L])/(zlim[2L]-zlim[1L])),0)])
#' 
#' @importFrom parallel detectCores makeCluster parSapply parLapply
#' 
NULL
#' 
#' @describeIn MCA
#' 
#' Main function to compute the multiscale codependence analysis
#' 
#' @export
MCA <- function(Y, X, emobj) {
  if (!inherits(emobj, "eigenmap"))
    stop("Parameter 'emobj' must be a 'eigenmap' object!")
  if (!is.matrix(Y))
    Y <- matrix(Y, length(Y), 1L, dimnames=list(names(Y), "Y"))
  if (!is.matrix(X))
    X <- matrix(X, length(X), 1L, dimnames=list(names(X), "X"))
  if (nrow(Y) != nrow(X))
    stop("Number of observations in Y and X do not match!")
  if (nrow(emobj$U) != nrow(Y))
    stop("Number of observations in Y does not match the number of lines in U.")
  Dnms <- list(Y=colnames(Y), X=colnames(X))
  mssd <- list(mY=NA, mX=NA, ssdY=NA, ssdX=NA)
  mssd[[1L]] <- colMeans(Y)
  mssd[[2L]] <- colMeans(X)
  YXc <- list(
    Yc = Y - rep(mssd$mY, each=nrow(Y)),
    Xc = X - rep(mssd$mX, each=nrow(X))
  )
  mssd[[3L]] <- colSums(YXc$Yc^2)
  mssd[[4L]] <- colSums(YXc$Xc^2)
  UpYXcb <- list(
    UpY = matrix(
      NA,
      ncol(emobj$U),
      ncol(Y),
      dimnames = list(colnames(emobj$U), Dnms$Y)
    ),
    UpX = matrix(
      NA,
      ncol(emobj$U),
      ncol(X),
      dimnames = list(colnames(emobj$U), Dnms$X)
    ),
    C = array(
      NA,
      dim = c(ncol(emobj$U),ncol(Y),ncol(X)),
      dimnames = list(colnames(emobj$U), Dnms$Y, Dnms$X)
    ),
    B = array(
      NA,
      dim=c(ncol(emobj$U),ncol(Y),ncol(X)),
      dimnames = list(colnames(emobj$U), Dnms$Y, Dnms$X)
    ),
    CM = matrix(
      NA,
      ncol(emobj$U),
      ncol(X),
      dimnames=list(colnames(emobj$U), Dnms$X)
    )
  )
  UpYXcb$UpY[] <- t(emobj$U) %*% YXc$Yc
  UpYXcb$UpX[] <- t(emobj$U) %*% YXc$Xc
  for (i in 1L:ncol(X)) {
    UpYXcb$C[,,i] <-
      UpYXcb$UpY/sqrt(rep(mssd[[3L]], each = ncol(emobj$U))) *
        UpYXcb$UpX[,i]/sqrt(mssd[[4L]][i])
    UpYXcb$B[,,i] <- UpYXcb$UpY/UpYXcb$UpX[,i]
    for (j in 1L:ncol(emobj$U)) {
      UpYXcb$CM[j,i] <-
        sqrt(sum((emobj$U[,j,drop=FALSE] %*%
                    UpYXcb$UpY[j,,drop = FALSE])^2) / sum(YXc$Yc^2)) *
        abs(UpYXcb$UpX[j,i] / sqrt(mssd[[4L]][i]))
    }
  }
  return(
    structure(
      list(
        data = list(Y = Y, X = X),
        emobj = emobj,
        UpYXcb = UpYXcb,
        test = NULL
      ),
      class = "cdp"
    )
  )
}
#' 
#' @describeIn MCA
#' 
#' Parametric statistical testing for multiscale codependence analysis
#' 
#' @export
test.cdp <- function(object, alpha = 0.05, max.step, response.tests = TRUE) {
  if (!inherits(object, "cdp"))
    stop("Parameter 'object' must be of class 'cdp'.")
  if (missing(max.step))
    max.step <- ncol(object$emobj$U)
  else
    max.step <- max.step[1L]
  us <- matrix(NA, nrow(object$emobj$U), 1L)
  uspY <- matrix(
    NA,
    1L,
    ncol(object$data$Y),
    dimnames = list(NULL, colnames(object$data$Y))
  )
  uspX <- matrix(
    NA,
    1L,
    ncol(object$data$X),
    dimnames = list(NULL, colnames(object$data$X))
  )
  Yc <- object$data$Y - rep(colMeans(object$data$Y), each=nrow(object$data$Y))
  Xc <- object$data$X - rep(colMeans(object$data$X), each=nrow(object$data$X))
  ord <- order(apply(object$UpYXcb$CM, 1L, max), decreasing=TRUE)
  bstX <- apply(object$UpYXcb$CM, 1L, which.max)
  ttable <- matrix(
      NA,
      0L,
      6L,
      dimnames = list(
        NULL,
        c("Variable","phi","df1","df2","Testwise p","Familywise p")
      )
    )
  if (response.tests) 
    respts <- array(
      numeric(0),
      dim = c(ncol(object$data$Y),4L,0L),
      dimnames = list(
        colnames(object$data$Y),
        c("tau","df","Testwise p","Familywise p"),
        NULL
      )
    )
  step <- 1L
  while (step != 0L) {
    us[] <- object$emobj$U[,ord[step]]
    uspY[] <- object$UpYXcb$UpY[ord[step],]
    uspX[] <- object$UpYXcb$UpX[ord[step],]
    df2 <- nrow(object$data$Y) - step - 1L
    Yhat <- us %*% uspY
    Xhat <- us %*% uspX
    Yc <- Yc - Yhat
    Xc <- Xc - Xhat
    phi_global <- df2^2 * sum(Yhat^2) *
      sum(Xhat[,bstX[ord[step]]]^2)/(sum(Yc^2) * sum(Xc[,bstX[ord[step]]]^2))
    ttable <-
      rbind(
        ttable,
        c(bstX[ord[step]],phi_global,ncol(object$data$Y),df2,NA,NA)
      )
    ttable[step,5L] <- pphi(
      phi_global,
      ncol(object$data$Y),
      df2,
      lower.tail = FALSE
    )
    ttable[step,6L] <-
      1 - (1 - ttable[step,5L])^((ncol(object$emobj$U) - step + 1) *
                                   ncol(object$data$X))
    if (response.tests) {
      tau_resp <- df2 * uspY[1L,] * uspX[bstX[ord[step]]] *
        (colSums(Yc^2) * sum(Xc[,bstX[ord[step]]]^2))^-0.5
      respts <-
        array(
          as.numeric(respts),
          dim = c(dim(respts)[1L],dim(respts)[2L],dim(respts)[3L] + 1L),
          dimnames = c(dimnames(respts)[1L:2],list(NULL))
        )
      respts[,1L,step] <- tau_resp
      respts[,2L,step] <- df2
      respts[,3L,step] <- 2 * ptau(abs(tau_resp), df2, lower.tail=FALSE)
      respts[,4L,step] <-
        1 - (1 - respts[,3L,step])^((ncol(object$emobj$U) - step + 1) *
                                      ncol(object$data$X))
    }
    if (ttable[step,6L] > alpha || step >= max.step) {
      rownames(ttable) <- colnames(object$emobj$U)[ord[1L:step]]
      if (response.tests)
        dimnames(respts)[[3L]] <- rownames(ttable)
      step <- 0
    }
    else step <- step + 1
  }
  signif <- list(U = ord[which(ttable[, 6L] <= alpha)])
  signif$X <- bstX[signif$U]
  return(
    structure(
      list(
        data = object$data,
        emobj = object$emobj,
        UpYXcb = object$UpYXcb,
        test = list(
          permute = FALSE,
          significant = signif,
          global = ttable,
          response = if (response.tests) respts else NULL,
          permutations = NULL)
        ),
      class = "cdp"
    )
  )
}
#' 
#' @describeIn MCA
#' 
#' Permutation testing for multiscale codependence analysis.
#' 
#' @export
permute.cdp <- function(object, permute, alpha = 0.05, max.step,
                        response.tests = TRUE) {
  if (!inherits(object, "cdp"))
    stop("Parameter 'object' must be of class 'cdp'.")
  if (missing(max.step))
    max.step <- ncol(object$emobj$U)
  else
    max.step <- max.step[1L]
  if (missing(permute))
    permute <- minpermute(
      alpha,
      ncol(object$emobj$U) * ncol(object$data$X),
      10L,
      3L
    )
  else
    permute <- permute[1L]
  us <- matrix(NA, nrow(object$emobj$U), 1L)
  uspY <- matrix(
    NA,
    1L,
    ncol(object$data$Y),
    dimnames = list(NULL, colnames(object$data$Y))
  )
  uspX <- matrix(
    NA,
    1L,
    ncol(object$data$X),
    dimnames = list(NULL, colnames(object$data$X))
  )
  Yc <- object$data$Y - rep(colMeans(object$data$Y), each=nrow(object$data$Y))
  Xc <- object$data$X - rep(colMeans(object$data$X), each=nrow(object$data$X))
  ord <- order(apply(object$UpYXcb$CM, 1L, max), decreasing=TRUE)
  bstX <- apply(object$UpYXcb$CM, 1L, which.max)
  ttable <-
    matrix(
      NA,
      0L,
      6L,
      dimnames = list(
        NULL,
        c("Variable","phi","df1","df2","Testwise p","Familywise p")
      )
    )
  perm_global <- matrix(
    NA,
    0L,
    2L,
    dimnames = list(NULL, c("phi*<phi","phi*>=phi"))
  )
  if (response.tests) {
    respts <- array(
      numeric(0),
      dim = c(ncol(object$data$Y),4L,0L),
      dimnames = list(
        colnames(object$data$Y),
        c("tau","df","Testwise p","Familywise p"),
        NULL
      )
    )
    perm_response <- array(
      numeric(0),
      dim=c(ncol(object$data$Y),3L,0L),
      dimnames = list(
        colnames(object$data$Y),
        c("tau*<=-|tau|","-|tau|<tau*<|tau|","tau*>=|tau|"),
        NULL
      )
    )
  }
  step <- 1L
  while (step != 0L) {
    us[] <- object$emobj$U[, ord[step]]
    uspY[] <- object$UpYXcb$UpY[ord[step],]
    uspX[] <- object$UpYXcb$UpX[ord[step],]
    df2 <- nrow(object$data$Y) - step - 1L
    Yhat <- us %*% uspY
    Xhat <- us %*% uspX
    Yc <- Yc - Yhat
    Xc <- Xc - Xhat
    phi_global0 <-
      sum(Yhat^2) * sum(Xhat[,bstX[ord[step]]]^2)/
      (sum(Yc^2) * sum(Xc[,bstX[ord[step]]]^2))
    ttable <- rbind(
      ttable,
      c(bstX[ord[step]],df2^2 * phi_global0,ncol(object$data$Y),df2,NA,NA)
    )
    perm_global <- rbind(perm_global, c(0L,1L))
    if (response.tests) {
      tau_resp0 <- uspY[1L,] * uspX[bstX[ord[step]]] *
        (colSums(Yc^2) * sum(Xc[, bstX[ord[step]]]^2))^-0.5
      respts <- array(
        as.numeric(respts),
        dim = c(dim(respts)[1L],dim(respts)[2L],dim(respts)[3L] + 1L),
        dimnames = c(dimnames(respts)[1L:2],list(NULL))
      )
      respts[,1L,step] <- df2 * tau_resp0
      respts[,2L,step] <- df2
      perm_response <-
        array(
          as.numeric(perm_response),
          dim=c(
            dim(perm_response)[1L],
            dim(perm_response)[2L],
            dim(perm_response)[3L] + 1L
          ),
          dimnames = c(dimnames(perm_response)[1L:2],list(NULL))
        )
      perm_response[,1L:2,step] <- 0
      perm_response[,3L,step] <- 1
      tmp <- .C(
        "mcapermute",
        as.double(phi_global0),
        as.double(abs(tau_resp0)),
        as.double(Yc),
        as.integer(ncol(Yc)),
        as.double(Xc[,bstX[ord[step]]]),
        as.double(us),
        as.integer(nrow(us)),
        perm_global = as.integer(perm_global[step,]),
        perm_response = as.integer(perm_response[,,step]),
        as.integer(permute),
        as.integer(TRUE)
      )
      perm_global[step,] <- tmp$perm_global
      perm_response[,,step] <- tmp$perm_response
    } else {
      perm_global[step,] <- .C(
        "mcapermute",
        as.double(phi_global0),
        as.double(),
        as.double(Yc),
        as.integer(ncol(Yc)),
        as.double(Xc[,bstX[ord[step]]]),
        as.double(us),
        as.integer(nrow(us)),
        perm_global = as.integer(perm_global[step,]),
        integer(),
        as.integer(permute),
        as.integer(FALSE)
      )$perm_global
    }
    ttable[step,5L] <- perm_global[step,2L] / (permute + 1)
    ttable[step,6L] <-
      1 - (1 - ttable[step, 5L])^((ncol(object$emobj$U) - step + 1) *
                                    ncol(object$data$X))
    if (response.tests) {
      respts[,3L,step] <-
        (perm_response[,1L,step] + perm_response[,3L,step]) / (permute + 1)
      respts[,4L,step] <-
        1 - (1 - respts[,3L,step])^((ncol(object$emobj$U) - step + 1) *
                                      ncol(object$data$X))
    }
    if (ttable[step, 6L] > alpha || step >= max.step) {
      rownames(ttable) <- colnames(object$emobj$U)[ord[1L:step]]
      if (response.tests) {
        dimnames(respts)[[3L]] <- rownames(ttable)
        dimnames(perm_response)[[3L]] <- rownames(ttable)
      }
      step <- 0
    }
    else step <- step + 1
  }
  signif <- list(U=ord[which(ttable[,6L] <= alpha)])
  signif$X <- bstX[signif$U]
  return(
    structure(
      list(
        data = object$data,
        emobj = object$emobj,
        UpYXcb = object$UpYXcb,
        test = list(
          permute = permute,
          significant = signif,
          global = ttable,
          response = if (response.tests) respts else NULL,
          permutations = list(
            global = perm_global,
            response = if (response.tests) perm_response else NULL
          )
        )
      ),
      class = "cdp"
    )
  )
}
#' 
#' @describeIn MCA
#' 
#' Permutation testing for multiscale codependence analysis using parallel
#' processing.
#' 
#' @export
parPermute.cdp <- function(object, permute, alpha = 0.05, max.step,
                           response.tests = TRUE, nnode, seeds,
                           verbose = TRUE, ...) {
  if (!inherits(object, "cdp"))
    stop("Parameter 'object' must be of class 'cdp'.")
  if (missing(max.step))
    max.step <- ncol(object$emobj$U)
  else
    max.step <- max.step[1L]
  if (missing(permute))
    permute <- minpermute(
      alpha,
      ncol(object$emobj$U) * ncol(object$data$X),
      10L,
      3L
    )
  else
    permute <- permute[1L]
  if (missing(nnode))
    nnode <- detectCores()
  if (verbose)
    cat("Starting a cluster of", nnode, "nodes... ")
  cl <- makeCluster(nnode, ...)
  nnode <- length(cl)
  if (verbose)
    cat("done.\n")
  if (nnode < 2L)
    warning(
      "Only a single worker could be recruited on that system.",
      "\nConsider using permute.mca() instead."
    )
  if (verbose)
    cat("Initializing workers:\n")
  if (missing(seeds))
    seeds <- as.integer(
      runif(nnode, -.Machine$integer.max, .Machine$integer.max)
    )
  if (verbose)
    cat("Random seeds given to workers:", seeds, "\n")
  parSapply(cl=cl, X=seeds, FUN=function(x) set.seed(x))
  parSapply(cl=cl, X=1L:nnode, function(x) require(codep))
  wpermute <- rep(permute%/%nnode, nnode)
  if (permute%%nnode)
    wpermute[1L:(permute%%nnode)] <- wpermute[1L:(permute%%nnode)]+1L
  C_mcapermute <- function(
    X,
    phi_global0,
    tau_ind0,
    rY,
    rx,
    us,
    perm_global,
    perm_response,
    ind
  ) .C(
    "mcapermute",
    as.double(phi_global0),
    as.double(abs(tau_ind0)),
    as.double(rY),
    as.integer(ncol(rY)),
    as.double(rx),
    as.double(us),
    as.integer(nrow(us)),
    as.integer(perm_global),
    as.integer(perm_response),
    as.integer(X),
    as.integer(ind)
  )[8L:9L]
  us <- matrix(NA, nrow(object$emobj$U), 1L)
  uspY <- matrix(
    NA,
    1L,
    ncol(object$data$Y),
    dimnames = list(NULL, colnames(object$data$Y))
  )
  uspX <- matrix(
    NA,
    1L,
    ncol(object$data$X),
    dimnames=list(NULL, colnames(object$data$X))
  )
  Yc <- object$data$Y - rep(colMeans(object$data$Y), each=nrow(object$data$Y))
  Xc <- object$data$X - rep(colMeans(object$data$X), each=nrow(object$data$X))
  ord <- order(apply(object$UpYXcb$CM, 1L, max), decreasing=TRUE)
  bstX <- apply(object$UpYXcb$CM, 1L, which.max)
  ttable <- matrix(
    NA,
    0L,
    6L,
    dimnames = list(
      NULL,
      c("Variable","phi","df1","df2","Testwise p","Familywise p")
    )
  )
  perm_global <- matrix(
    NA,
    0L,
    2L,
    dimnames = list(NULL, c("phi*<phi","phi*>=phi"))
  )
  if (response.tests) {
    respts <- array(
      numeric(0),
      dim = c(ncol(object$data$Y),4L,0L),
      dimnames = list(
        colnames(object$data$Y),
        c("tau","df","Testwise p","Familywise p"),
        NULL
      )
    )
    perm_response <- array(
      numeric(0),
      dim=c(ncol(object$data$Y),3L,0L),
      dimnames=list(
        colnames(object$data$Y),
        c("tau*<=-|tau|","-|tau|<tau*<|tau|","tau*>=|tau|"),
        NULL
      )
    )
  }
  if (verbose)
    cat("Performing permutation tests")
  step <- 1L
  while (step != 0L) {
    us[] <- object$emobj$U[,ord[step]]
    uspY[] <- object$UpYXcb$UpY[ord[step],]
    uspX[] <- object$UpYXcb$UpX[ord[step],]
    df2 <- nrow(object$data$Y) - step - 1L
    Yhat <- us %*% uspY
    Xhat <- us %*% uspX
    Yc <- Yc - Yhat
    Xc <- Xc - Xhat
    phi_global0 <- sum(Yhat^2) * sum(Xhat[,bstX[ord[step]]]^2)/
      (sum(Yc^2) * sum(Xc[,bstX[ord[step]]]^2))
    ttable <- rbind(
      ttable,
      c(bstX[ord[step]],df2^2 * phi_global0, ncol(object$data$Y),df2,NA,NA)
    )
    perm_global <- rbind(perm_global, c(0L,0L))
    if (response.tests) {
      tau_resp0 <- uspY[1L,] * uspX[bstX[ord[step]]] *
        (colSums(Yc^2) * sum(Xc[,bstX[ord[step]]]^2))^-0.5
      respts <- array(
        as.numeric(respts),
        dim = c(dim(respts)[1L],dim(respts)[2L],dim(respts)[3L] + 1L),
        dimnames = c(dimnames(respts)[1L:2L],list(NULL))
      )
      respts[,1L,step] <- df2 * tau_resp0
      respts[,2L,step] <- df2
      perm_response <- array(
        as.numeric(perm_response),
        dim = c(
          dim(perm_response)[1L],
          dim(perm_response)[2L],
          dim(perm_response)[3L] + 1L
        ),
        dimnames = c(dimnames(perm_response)[1L:2L],list(NULL))
      )
      perm_response[,,step] <- 0
      wtmp <- parLapply(
        cl = cl,
        X = wpermute,
        fun = C_mcapermute,
        phi_global0 = phi_global0,
        tau_ind0 = tau_resp0,
        rY = Yc,
        rx = Xc[,bstX[ord[step]]],
        us = us,
        perm_global = perm_global[step,],
        perm_response = perm_response[,,step],
        ind = TRUE
      )
      tmp <- wtmp[[1L]]
      if (nnode > 1L)
        for (i in 2L:nnode) {
          tmp[[1L]] <- tmp[[1L]] + wtmp[[i]][[1L]]
          tmp[[2L]] <- tmp[[2L]] + wtmp[[i]][[2L]]
        }
      perm_global[step,] <- tmp[[1L]]
      perm_global[step,2L] <- perm_global[step,2L] + 1
      perm_response[,,step] <- tmp[[2L]]
      perm_response[,3L,step] <- perm_response[,3L,step] + 1
    }
    else {
      wtmp <- parLapply(
        cl = cl,
        X = wpermute,
        fun = C_mcapermute,
        phi_global0 = phi_global0,
        tau_ind0 = tau_resp0,
        rY = Yc,
        rx = Xc[,bstX[ord[step]]],
        us = us,
        perm_global = perm_global[step,],
        perm_response = integer(),
        ind = FALSE
      )
      tmp <- wtmp[[1L]][[1L]]
      if(nnode > 1L)
        for(i in 2L:nnode)
          tmp <- tmp + wtmp[[i]][[1L]]
      perm_global[step,] <- tmp
      perm_global[step,2L] <- perm_global[step,2L] + 1
    }
    ttable[step,5L] <- perm_global[step,2L]/(permute + 1)
    ttable[step,6L] <-
      1 - (1 - ttable[step,5L])^((ncol(object$emobj$U) - step + 1)*
                                   ncol(object$data$X))
    if (response.tests) {
      respts[,3L,step] <-
        (perm_response[,1L,step] + perm_response[,3L,step])/(permute + 1)
      respts[,4L,step] <-
        1 - (1 - respts[,3L,step])^((ncol(object$emobj$U) - step + 1)*
                                      ncol(object$data$X))
    }
    if (verbose)
      cat(".")
    if (ttable[step,6L] > alpha || step >= max.step) {
      rownames(ttable) <- colnames(object$emobj$U)[ord[1L:step]]
      if (response.tests) {
        dimnames(respts)[[3L]] <- rownames(ttable)
        dimnames(perm_response)[[3L]] <- rownames(ttable)
      }
      step <- 0
    } else step <- step + 1
  }
  if (verbose)
    cat("done.\nStopping cluster.\n")
  parallel::stopCluster(cl)
  signif <- list(U = ord[which(ttable[, 6L] <= alpha)])
  signif$X <- bstX[signif$U]
  return(
    structure(
      list(
        data = object$data,
        emobj = object$emobj,
        UpYXcb = object$UpYXcb,
        test = list(
          permute = permute,
          significant = signif,
          global = ttable,
          response = if (response.tests) respts else NULL,
          permutations = list(
            global = perm_global,
            response = if (response.tests) perm_response else NULL
          )
        ),
        seeds = seeds
      ),
      class = "cdp"
    )
  )
}
#' 
