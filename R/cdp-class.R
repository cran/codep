## **************************************************************************
##
##    (c) 2018-2023 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    **cdp-class definition**
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
#' 
#' Class and Methods for Multiscale Codependence Analysis (MCA)
#' 
#' A class and set of methods to handle the results of Multiscale Codependence
#' Analysis.
#' 
#' @name cdp-class
#' 
#' @docType class
#' 
#' @param x A \code{\link{cdp-class}} object.
#' @param object A \code{\link{cdp-class}} object.
#' @param col A vector of color values to be used for plotting the multivariate
#' codependence coefficients.
#' @param col.signif Color of the frame used to mark the statistically
#' significant codependence coefficients.
#' @param main Text for the main title of the plot.
#' @param selection A numeric vector of indices or character vector variable
#' names to test or force-use. Mandatory if \code{object} is untested.
#' @param components A boolean specifying whether the components of fitted or
#' predicted values associated with single eigenfunctions in the map should be
#' returned.
#' @param newdata A list with elements $X, $meanY, and $target that contain the
#' information needed to make predictions (see details).
#' @param ... Further parameters to be passed to other functions or methods.
#' 
#' @details
#' The `fitted`, `residuals`, and `predict` methods return a matrix of fitted,
#' residuals, or predicted values, respectively. The `fitted` and `predict`
#' methods return a list a list when argument `component` is \code{TRUE}. The
#' list contains the `fitted` or `predicted` values as a first element and an
#' array `components` as a second. That 3-dimensional array has one matrix for
#' each statistically significant codependence coefficient.
#' 
#' For making predictions, argument \code{newdata} may contain three elements:
#' `$X`, a matrix of new values of the explanatory variables, `$meanY`, a vector
#' of the predicted mean values of the responses, and `$target`, a matrix of
#' target scores for arbitrary locations within the study area. When no `$X` is
#' supplied, the descriptor given to \code{\link{MCA}} is recycled, while when
#' no `$meanY` is supplied, the mean values of the response variables given to
#' \code{\link{MCA}} are used.
#' 
#' Finally, when element `$target` is omitted from argument \code{newdata},
#' predictions are made at the sites were observations were done. When none of
#' the above is provided, or if \code{newdata} is omitted when calling the
#' prediction method, the behaviour of the `predict` method is identical to
#' that of the `fitted` method.
#' 
#' From version 0.7-1, \link{cdp-class} replaces the former class
#' \code{mca} used by \link{codep-package} because the standard package MASS
#' also had S3 methods for a class named \code{mca} that were overwritten by
#' those of \link{codep-package}.
#' 
#' @format cdp-class objects contain:
#' \describe{
#'   \item{data}{ A list with two elements: the first being a copy of the
#'     response (`Y`) and the second being a copy of the explanatory variables
#'     (`X`). This is the variables that were given to \code{\link{MCA}}. }
#'   \item{emobj}{ The \link{eigenmap-class} object that was given to
#'     \code{\link{MCA}}. }
#'   \item{UpYXcb}{ A list with five elements: the first (`UpY`) is a matrix of
#'     the cross-products of structuring variable (`U`) and the response
#'     variable `Y`, the second (`UpX`) is a matrix of the cross-product of the
#'     structuring variable and the explanatory variables (`X`), the third (`C`)
#'     is a 3-dimensional array of the codependence coefficients, the fourth
#'     (`B`) is a 3-dimensional array of the coregression coefficients, and the
#'     fifth (`CM`) is a matrix of the multivariate codependence coefficients. }
#'   \item{test}{ Results of statistical testing as performed by
#'     \code{\link{test.cdp}} or \code{\link{permute.cdp}}. \code{NULL} if no
#'     testing was performed, such as when only \code{\link{MCA}} had been
#'     called. The results of statistical testing is a list containing the
#'     following members:
#'     \describe{
#'       \item{$permute}{ The number of randomized permutations used by
#'         \code{permute.cdp} for permutation testing. 0 or \code{FALSE} for
#'         parametric testing obtained using \code{\link{test.cdp}}. }
#'       \item{$significant}{ The indices of codependence coefficient describing
#'         statistically significant codependence between `Y` and `X`, in
#'         decreasing order of magnitude. }
#'       \item{$global}{ The testing table (a 5-column matrix) with phi
#'         statistics, degrees-of-freedom, and testwise and familywise
#'         probabilities of type I (alpha) error. It contains one line for each
#'         statistically significant global coefficient (if any) in addition to
#'         test results for the first, non-significant coefficient, on which the
#'         testing procedure stopped. }
#'       \item{$response}{ Tests of every single response variable (a
#'         3-dimensional array), had such tests been requested while calling the
#'         testing function, \code{NULL} otherwise. }
#'       \item{$permutations}{ Details about permutation testing not shown in
#'         `test$global` or `test$response`. \code{NULL} for parametric
#'         testing. }
#'     }
#'   }
#' }
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
#' @importFrom graphics axis box hist image par rect
#' @importFrom stats runif
#' 
NULL
#' 
#' @describeIn cdp-class
#' 
#' Print method for cdp-class objects.
#' 
#' @export
print.cdp <- function (x, ...) {
  cat(
    "\nMultiple Multi-scale Codependence Analysis\n",
    "---------------------------\n\n"
  )
  cat(
    ncol(x$data$X),
    " explanatory variable",
    if (ncol(x$data$X) > 1L) "s",
    "\n\n",
    sep = ""
  )
  print(
    signif(
      cbind(
        x$emobj$lambda,
        x$UpYXcb$CM
      ),
      4
    )
  )
  if (!is.null(x$test)) {
    cat("--------------------\nGlobal testing information is available\n")
    if (!is.null(x$test$response)) 
      cat(
        "Hypothesis test",
        if (ncol(x$data$X)>1L) "s",
        " also available for the response",
        if (ncol(x$data$X) > 1L) "s",
        "\n",
        sep = ""
      )
    else
      cat(
        "Individual test",
        if (ncol(x$data$X) > 1L) "s",
        " unavailable\n",
        sep = ""
      )
  } else
    cat("\n")
  return(invisible(NULL))
}
#' 
#' @describeIn cdp-class
#' 
#' Plot method for cdp-class objects.
#' 
#' @export
plot.cdp <- function (x, col, col.signif = 2, main = "", ...) {
  if(missing(col))
    col <- grey(seq(1, 0, length.out=256))
  mar <- par()$mar
  z <- log10(x$UpYXcb$CM + 1e-04)
  par(
    mar = c(mar[1L],mar[2L],mar[3L],0.75),
    fig = c(0,0.875 - 0.025 * (mar[4L] - 2.1),0,1)
  )
  image(
    y = 1L:ncol(x$data$X),
    x = 1L:ncol(x$emobj$U),
    z = z,
    zlim = c(-4, 1e-04),
    col = col,
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = main,
    ...
  )
  box(...)
  axis(1, at=1L:ncol(x$emobj$U), labels=colnames(x$emobj$U), ...)
  axis(2, at=1L:ncol(x$data$X), labels=colnames(x$data$X), ...)
  if (!is.null(x$test$signif)) 
    rect(
      xleft = x$test$signif$U - 0.5,
      xright = x$test$signif$U + 0.5,
      ybottom = x$test$signif$X - 0.5,
      ytop = x$test$signif$X + 0.5,
      border = col.signif,
      density = NULL,
      ...
    )
  par(
    mar = c(mar[1L],0.75,mar[3L],mar[4L]),
    fig = c(0.875 - 0.025 * (mar[4L] - 2.1),1,0,1),
    new = TRUE
  )
  image(
    z = matrix(seq(-4, 1e-04, length.out=256L), 1L, 256L),
    x = 0,
    y = seq(-4, 1e-04, length.out=256L),
    col = col,
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = "",
    ...
  )
  box(...)
  axis(4, labels=10^seq(-4, 0, 1), at=seq(-4, 0, 1), ...)
  par(mar=mar, fig=c(0,1,0,1))
  return(invisible(NULL))
}
#' 
#' @describeIn cdp-class
#' 
#' Summary method for cdp-class objects.
#' 
#' @export
summary.cdp <- function(object, ...) {
  cat(
    "\nMultiple Multi-scale Codependence Analysis\n",
    "---------------------------\n\n"
  )
  cat(
    ncol(object$data$X),
    " explanatory variable",
    if (ncol(object$data$X)>1) "s",
    "\n\n",
    sep=""
  )
  if (is.null(object$test)) {
    cat("\nNo testing informations available\n\n")
  } else {
    cat("\nTest table:\n")
    tmp <- data.frame(
      object$test$global[,1L:4,drop=FALSE],
      character(nrow(object$test$global)),
      stringsAsFactors = FALSE
    )
    colnames(tmp)[5L] <- "Pr(>|phi|)"
    tmp[,"Variable"] <- colnames(object$data$X)[tmp[,"Variable"]]
    tmp[object$test$global[,6L] < 2.2e-16,"Pr(>|phi|)"] <- "<2.2e-16"
    tmp[object$test$global[,6L] >= 2.2e-16,"Pr(>|phi|)"] <-
      signif(object$test$global[object$test$global[,6L] >= 2.2e-16,6L], 2L)
    print(tmp)
    cat("\n")
    if (is.null(object$test$response))
      cat("No individual response testing available\n")
    else {
      cat("Individual response tests:\n")
      tmp <- matrix(
        NA,
        nrow(object$test$global),
        6L,
        dimnames = list(
          rownames(object$test$global),
          NULL
        )
      )
      for (i in 1L:nrow(object$test$global))
        tmp[i,] <- hist(
          x = object$test$response[, 4L, i],
          plot = FALSE,
          breaks = c(0,1e-04,0.001,0.01,0.05,0.1,1)
        )$counts
      tmp <- data.frame(
        tmp,
        Variable = colnames(object$data$X)[object$test$global[,1L]]
      )
      colnames(tmp) <- c(
        "p<=0.0001",
        "0.0001>p>=0.001",
        "0.001>p>=0.01",
        "0.01>p>=0.05",
        "0.05>p>=0.1",
        "p>0.1",
        "Variable"
      )
      print(tmp)
      cat("\n")
    }
  }
  return(invisible(TRUE))
}
#' 
#' @describeIn cdp-class
#' 
#' Fitted method for cdp-class objects.
#' 
#' @export
fitted.cdp <- function (object, selection, components = FALSE, ...) {
  if(missing(selection)) {
    if (!is.null(object$test))
      selection <- object$test$significant
    else
      stop(
        "No testing informations available: user must identify ",
        "the relevant coefficients."
      )
  } else {
    if (is.null(selection$U))
      stop("Parameter 'selection' must be a list with an element $U")
  }
  if (components) {
    cpns <- array(
      NA,
      dim = c(nrow(object$data$Y),ncol(object$data$Y),length(selection$U))
    )
    dimnames(cpns) <- list(
      rownames(object$data$Y),
      colnames(object$data$Y),
      colnames(object$emobj$U)[selection$U]
    )
  }
  if (length(selection$U)) {
    by <- object$UpYXcb$UpY[selection$U,,drop=FALSE]
    fit <- object$emobj$U[,selection$U,drop=FALSE] %*% by +
      rep(colMeans(object$data$Y), each = nrow(object$data$Y))
    if (components) 
      for (i in 1L:length(selection$U)) {
        cpns[,,i] <- object$emobj$U[,selection$U[i],drop = FALSE] %*%
          by[i,,drop = FALSE]
      }
  }
  if (components)
    return(list(fitted=fit, components=cpns))
  else
    return(fit)
}
#' 
#' @describeIn cdp-class
#' 
#' Residuals method for cdp-class objects.
#' 
#' @export
residuals.cdp <- function (object, selection, ...) {
  if(missing(selection)) {
    if (!is.null(object$test))
      selection <- object$test$significant
    else
      stop(
        "No testing informations available: user must identify ",
        "the relevant coefficients."
      )
  } else {
    if (is.null(selection$U))
      stop("Parameter 'selection' must be a list with an element $U")
  }
  res <- object$data$Y - rep(colMeans(object$data$Y), each=nrow(object$data$Y))
  if (length(selection$U)) {
    by <- object$UpYXcb$UpY[selection$U,]
    res <- res - object$emobj$U[,selection$U,drop=FALSE] %*% by
  }
  return(res)
}
#' 
#' @describeIn cdp-class
#' 
#' Predict method for cdp-class objects.
#' 
#' @export
predict.cdp <- function (object, selection, newdata, components = FALSE, ...) {
  if (missing(newdata))
    return(fitted.cdp(object, selection = selection))
  if(missing(selection)) {
    if (!is.null(object$test))
      selection <- object$test$significant
    else
      stop(
        "No testing informations available: user must identify ",
        "the relevant coefficients."
      )
  } else {
    if (is.null(selection$U)||is.null(selection$X))
      stop("Parameter 'selection' must be a list with elements $U and $X")
  }
  if(!is.null(newdata$X)) {
    if ((NROW(newdata$X) != nrow(object$data$X)) ||
        (NCOL(newdata$X) != ncol(object$data$X)))
      stop(
        "'newdata$X' (",
        NROW(newdata$X),
        "x",
        NCOL(newdata$X),
        ") is not compatible with the original descriptors",
        nrow(object$data$X),
        "x",
        ncol(object$data$X)
      )
    if (NROW(newdata$X) != nrow(object$emobj$U))
      stop(
        "The number of observations in 'newdata$X' does not match ",
        "the number of observations."
      )
    by <- matrix(
      NA,
      length(selection$U),
      ncol(object$data$Y),
      dimnames = list(
        colnames(object$emobj$U)[selection$U],
        colnames(object$data$Y)
      )
    )
    for (i in 1L:length(selection$X)) {
      by[i,] <- object$UpYXcb$B[selection$U[i],,selection$X[i]] *
        as.numeric(
          t(object$emobj$U[,selection$U[i],drop=FALSE]) %*%
            (newdata$X[,selection$X[i]] - mean(newdata$X[,selection$X[i]]))
        )
    }
  } else
    by <- object$UpYXcb$UpY[selection$U,,drop=FALSE]
  if (is.null(newdata$meanY))
    newdata$meanY <- colMeans(object$data$Y)
  else
    if (length(newdata$meanY) != ncol(object$data$Y))
      stop(
        "The number of means in 'newdata$meanY' does not match ",
        "the number of response variable."
      )
  if (is.null(newdata$target))
    newdata$target <- object$emobj$U
  else
    if (ncol(newdata$target) != ncol(object$emobj$U))
      stop(
        "Incorrect number of eigenfunctions (columns) in 'newdata$meanY': ",
        ncol(newdata$target),
        ", while ",
        ncol(object$emobj$U),
        " is expected"
      )
  if (components) {
    cpns <- array(
      NA,
      dim = c(nrow(newdata$target),ncol(object$data$Y),length(selection$U))
    )
    dimnames(cpns) <- list(
      rownames(newdata$target),
      colnames(object$data$Y),
      colnames(newdata$target)[selection$U]
    )
  }
  pred <- newdata$target[,selection$U,drop=FALSE] %*% by +
    rep(newdata$meanY, each=nrow(newdata$target))
  if (components) {
    for (i in 1L:length(selection$U)) {
      cpns[,,i] <- newdata$target[,selection$U[i],drop=FALSE] %*%
        by[i,,drop=FALSE]
    }
    return(list(predicted=pred, components=cpns))
  } else return(pred)
}
#' 
