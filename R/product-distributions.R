## **************************************************************************
##
##    (c) 2018-2023 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    **Product distribution functions**
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
#' Frequency Distributions for MCA Parametric Testing
#' 
#' Density and distribution functions of the phi statistic, which is the product
#' of two Fisher-Snedecor distributions or the tau statistic, which is the
#' product of two Student's t distributions.
#' 
#' @name product-distribution
#' 
#' @param x,q A vector of quantile.
#' @param nu1,nu2,nu Degrees of freedom (>0, may be non-integer; \code{Inf} is
#' allowed.
#' @param lower.tail Logical; if TRUE (default), probabilities are P[X <= x],
#' otherwise, P[X > x].
#' @param tol The tolerance used during numerical estimation.
#' 
#' @return \code{dphi} and \code{dtau} return the density functions, whereas
#' \code{pphi} and \code{ptau} return the distribution functions.
#' 
#' @details The density distribution of a variable \code{z} that is the product
#' of two random variables `x` and `y` with density distributions f(x) and g(y),
#' respectively, is the integral f(x) * g(z/x) / abs(x) dx over the intersection
#' of the domains of `x` and `y`.
#' 
#' \code{dphi} estimates density values using numerical integration
#' (\code{\link{integrate}}) the Fisher-Scedecor \code{\link{df}} density
#' distribution function. Following the algebra of Multiscale Codependence
#' Analysis, f(x) has df1 = nu1 and df2 = nu1 * nu2 degrees of freedom and g(x)
#' has 'df1 = 1' and 'df2 = nu2' degrees of freedom. Hence, that product
#' distribution has two parameters.
#' 
#' \code{pphi} integrates \code{dphi} in the interval [0,q] when
#' \code{lower.tail = TRUE} (the default) and on the interval [q,\code{Inf}]
#' when \code{lower.tail = FALSE}.
#' 
#' \code{dtau} and \code{ptau} are similar to \code{dphi} and \code{pphi},
#' respectively. \code{pphi} integrates \code{dphi}, with f(x) and f(y) being
#' two Student's t distribution with \code{nu} degrees of freedom. It is called
#' by functions \code{\link{test.cdp}} and \code{\link{permute.cdp}} to perform
#' hypothesis tests for single response variables, in which case unilateral
#' tests can be performed.
#' 
#' @author \packageAuthor{codep}
#' Maintainer: \packageMaintainer{codep}
#' 
#' @seealso \link{test.cdp}
#' 
#' @references
#' Springer, M. D. 1979. The algebra of random variables. John Wiley and Sons
#' Inc., Hoboken, NJ, USA.
#' 
#' Guénard, G., Legendre, P., Boisclair, D., and Bilodeau, M. 2010. Multiscale
#' codependence analysis: an integrated approach to analyse relationships across
#' scales. Ecology 91: 2952-2964
#' 
#' Guénard, G. Legendre, P. 2018. Bringing multivariate support to multiscale
#' codependence analysis: Assessing the drivers of community structure across
#' spatial scales. Meth. Ecol. Evol. 9: 292-304
#' 
#' @examples
#' ### Displays the phi probability distribution for five different numbers
#' ### of degrees of freedom:
#' x <- 10^seq(-4, 0.5, 0.05)
#' plot(y = dphi(x, 1, 10), x = x, type = "l", col = "black", las = 1,
#' ylab = "pdf", ylim = c(0, 0.5))
#' lines(y = dphi(x, 3, 10), x = x, col = "purple")
#' lines(y = dphi(x, 5, 70), x = x, col = "blue")
#' lines(y = dphi(x, 12, 23), x = x, col = "green")
#' lines(y = dphi(x, 35, 140), x = x, col = "red")
#' 
#' ### Displays the density distribution function for 10 degrees of freedom.
#' x <- 10^seq(-4, 0.5, 0.05)
#' y <- dphi(x, 5, 70)
#' plot(y = y, x = x, type = "l", col = "black", las = 1, ylab = "Density",
#'      ylim = c(0, 0.5))
#' polygon(x = c(x[81L:91], x[length(x)], 1), y = c(y[81L:91], 0, 0),
#'         col = "grey")
#' text(round(pphi(1, 5, 70, lower.tail=FALSE), 3), x = 1.75, y = 0.05)
#' 
#' ## Idem for the tau distribution:
#' x <- c(-(10^seq(0.5, -4, -0.05)), 10^seq(-4, 0.5, 0.05))
#' plot(y = dtau(x, 1), x = x, type = "l", col = "black", las = 1,
#'      ylab = "pdf", ylim = c(0, 0.5))
#' lines(y = dtau(x, 2), x = x, col = "purple")
#' lines(y = dtau(x, 5), x = x, col="blue")
#' lines(y = dtau(x, 10), x = x, col="green")
#' lines(y = dtau(x, 100), x = x, col="red")
#' 
#' y <- dtau(x, 10)
#' plot(y = y, x = x, type = "l", col = "black", las = 1, ylab = "Density",
#'      ylim = c(0, 0.5))
#' polygon(x = c(x[which(x==1):length(x)], x[length(x)],1),
#'         y = c(y[which(x==1):length(x)], 0, 0), col = "grey")
#' text(round(ptau(1, 10, lower.tail = FALSE), 3), x = 1.5, y = 0.03)
#' polygon(x = c(-1, x[1L], x[1L:which(x==-1)]),
#'         y = c(0, 0, y[1L:which(x==-1)]), col="grey")
#' text(round(ptau(-1, 10), 3), x = -1.5, y = 0.03)
#' 
#' @importFrom stats df dt integrate
#' 
NULL
#' 
#' @describeIn product-distribution
#' 
#' Probability density function for the phi statistics
#' 
#' @export
dphi <- function(x, nu1, nu2, tol = .Machine$double.eps ^ 0.5) {
  res <- numeric(length(x))
  nu1 <- rep(nu1, length.out=length(x))
  nu2 <- rep(nu2, length.out=length(x))
  ## Integrand function to be integrated over positive real numbers.
  f <- function(z, x, nu1, nu2)
    df(z, nu1, nu1 * nu2) * df(x / z, 1, nu2) / z
  ## abs(z) always == z for real positive.
  ## The integrand is strictly positive.
  for (i in 1L:length(x))
    res[i] <- integrate(
      f,
      lower = 0,
      upper = Inf,
      x = x[i],
      nu1 = nu1[i],
      nu2 = nu2[i],
      rel.tol = tol
    )$value
  return(res)
}
#'
#' @describeIn product-distribution
#' 
#' Distribution function for the phi statistics
#' 
#' @export
pphi <- function(q, nu1, nu2, lower.tail = TRUE,
                 tol = .Machine$double.eps^0.5) {
  res <- numeric(length(q))
  nu1 <- rep(nu1, length.out=length(q))
  nu2 <- rep(nu2, length.out=length(q))
  if (lower.tail)
    for (i in 1L:length(q))
      res[i] <- 1 - integrate(
        dphi,
        lower = q[i],
        upper = Inf,
        nu1 = nu1[i],
        nu2 = nu2[i],
        tol = tol,
        rel.tol = tol
      )$value
  else
    for (i in 1L:length(q))
      res[i] <- integrate(
        dphi,
        lower = q[i],
        upper = Inf,
        nu1 = nu1[i],
        nu2 = nu2[i],
        tol = tol,
        rel.tol = tol
      )$value
  return(res)
}
#'
#' @describeIn product-distribution
#' 
#' Probability density function for the tau statistics
#' 
#' @export
dtau <- function(x, nu, tol = .Machine$double.eps ^ 0.5) {
  res <- numeric(length(x))
  nu <- rep(nu, length.out=length(x))
  ## Integrand function to be integrated over positive real numbers.
  f <- function(z, x, nu)
    dt(z, nu) * dt(x / z, nu) / z   ## abs(z) always == z for real positive.
  ## The integrand being symmetric, only the positive is integrated and the
  ## result multiplied by two.
  for (i in 1L:length(x))
    res[i] <- 2 * integrate(
      f,
      lower = 0,
      upper = Inf,
      x = x[i],
      nu = nu[i],
      rel.tol = tol
    )$value
  return(res)
}
#'
#' @describeIn product-distribution
#' 
#' Distribution function for the tau statistics
#' 
#' @export
ptau <- function(q, nu, lower.tail = TRUE, tol = .Machine$double.eps^0.5) {
  res <- rep(0.5, length(q))
  nu <- rep(nu, length.out=length(q))
  ## Only the positive 
  signq <- sign(q)
  q <- abs(q)
  ## Code avoid calculate PDF(0) because the function has a singularity at that
  ## point.
  if (lower.tail) {
    for (i in which(!!signq)) {
      if (signq[i]==1)
        res[i] <- 1 - integrate(
          dtau,
          lower = q[i],
          upper = Inf,
          nu = nu[i],
          tol = tol,
          rel.tol = tol
        )$value
      else
        res[i] <- integrate(
          dtau,
          lower = q[i],
          upper = Inf,
          nu = nu[i],
          tol = tol,
          rel.tol = tol
        )$value
    }
  } else {
    for (i in which(!!signq)) {
      if (signq[i] == 1)
        res[i] <- integrate(
          dtau,
          lower = q[i],
          upper = Inf,
          nu = nu[i],
          tol = tol,
          rel.tol = tol
        )$value
      else
        res[i] <- 1 - integrate(
          dtau,
          lower = q[i],
          upper = Inf,
          nu = nu[i],
          tol = tol,
          rel.tol = tol
        )$value
    }
  }
  return(res)
}
#' 
