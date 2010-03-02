% --- Source file: ./man/mevmap.Rd ---
% File src/library/codep/man/mevmap.Rd
% Part of the R package, http://www.R-project.org
% Copyright 2009 Guenard, Guillaume and Legendre, Pierre
% Distributed under GPL 2 or later
\name{mca}
\alias{mca}
\alias{test.mca}
\alias{permute.mca}
\alias{print.mca}
\alias{summary.mca}
\alias{plot.mca}
\alias{fitted.mca}
\alias{residuals.mca}
\alias{predict.mca}
\title{ Multiscale Codependence Analysis }
\description{ Functions and methods to handle Multiscale Codependence
  Analysis. }
\usage{
mca(y, x, memobj)
test.mca(mcaobj, alpha = 0.05, max.step = Inf)
permute.mca(mcaobj, permute = NA, alpha = 0.05, max.step = Inf)
\method{print}{mca}(x, ...)
\method{summary}{mca}(object, ...)
\method{plot}{mca}(x, ...)
\method{fitted}{mca}(object, which=NA, components=FALSE, ...)
\method{residuals}{mca}(object, which=NA, ...)
\method{predict}{mca}(object, which=NA, newdata=NA, components=FALSE, ...)
}
\arguments{
  \item{y}{ A numeric vector containing the explained variable. }
  \item{x}{ A numeric vector containing the explanatory variable for
    functions \code{test.mca} or \code{permute.mca} or a mca object for
    methods \code{print} and \code{plot}. }
  \item{memobj}{ An object of \code{\link{class}} \sQuote{mem} obtained
    using \code{\link{meigmap}} or \code{\link{mem}}. }
  \item{object, mcaobj}{ An object of \code{\link{class}} \sQuote{mca}
    obtained from \code{mca}, \code{test.mca}, or \code{permute.mca} }
  \item{alpha}{ The type I (\eqn{\alpha}{alpha}) error threshold used by
    the testing procedure. }
  \item{max.step}{ The maximum number of steps to perform when testing
    for statistical significance. }
  \item{permute}{ The number of random permutations used for
    testing. Obtained through \code{\link{minpermute}} when let to
    \code{NA}. }
  \item{which}{ A numeric vector of indices or character vector variable
    names to force-use. Mandatory if \code{mcaobj} is untested. }
  \item{components}{ A boolean specifying whether the components of
    fitted or predicted values associated with single eigenvector in the
    map (see \code{\link{mem}}) should be returned. }
  \item{newdata}{ A numeric vector containing new values of the
    explonatory variable. }
  \item{...}{ Further parameters to be passed to other functions or
    methods (currently ignored). }
}
\details{
  Multiscale Codependence Analysis (MCA) allows to calculate
  correlation-like (i.e.codependence) coefficients between two variables
  with respect to structuring variables (Moran's eigenvector maps). The
  analysis is performed using function mca. The purpose of this function
  is limited to parameter fitting. Test procedures are handled through
  test.mca (parametric testing) or permute.mca (permutation
  testing). Additionaly, methods are provided for printing, obtaining
  testing summary, plotting results, calculating fitted and residuals
  values, and making predictions. }
\value{
  An object of \code{\link{class}} \sQuote{mca} is a list containing:
  \item{data}{ A copy of the explaned (\code{y}) and explanatory
    (\code{x}) variables that were given to \code{mca}. }
  \item{memobj}{ The object of \code{\link{class}} \sQuote{mem} that was
    given to \code{mca}. }
  \item{Upyxcb}{ A 4 columns matrix containing the vectors of
    cross-products of structuring variable (\code{U} see
    \code{\link{mem}}) with variable \code{y} and \code{x}, the
    codependence coefficients \eqn{C}{C} and the coregression coefficients
    \eqn{B}{B}. }
  \item{test}{ Results of statistical testing. \code{NULL} if no testing
    was performed. The results of statistical testing is a list
    containing: }  
  \item{permute}{ The number of randomized permutations used by
    \code{permute.mca} for permutation testing. 0 or \code{FALSE} for
    parametric testing obtained using \code{test.mca}. }
  \item{significant}{ The indices of codependence coefficient describing
    statistically significant codependence between \code{y} and
    \code{x}, in decreasing order of importance. }
  \item{test.table}{ The testing table (a 4 columns matrix) with tau
    statistics, degrees-of-freedom, and testwise and familywise
    probabilities of type I (\eqn{\alpha}{alpha})) error. It contains one line for
    each statistically significant coefficient (if any) and an
    additional line showing the last tested (i.e. non-significant)
    coefficient. }
  \item{details}{ Details about permutation testing not shown in
    \code{test.table}. \code{NULL} for parametric testing. }
  The \code{fitted}, \code{residuals}, and \code{predict} methods return
  a single-column matrix of fitted, residuals, or predicted values,
  respectively. The \code{fitted} and \code{predict} methods return a
  list a list when the parameter \code{component} is \code{TRUE}. The
  list contains the \code{fitted} or \code{predicted} values as a first
  item and a matrix \code{components} as a second. This matrix has one
  column for each statistically significant codependence coefficient.
}
%\references{  }
\author{ Guillaume Guenard, Laboratoire evolution et diversite
  biologique, CNRS / Universite Paul-Sabatier, Toulouse, France. }
\seealso{
  \code{\link{meigmap}}
}
\examples{
#
### Exemple 1: Atlantic salmon in the St-Marguerite river, Quebec, Canada.
#
data(Salmon)
map <- meigmap(x=Salmon[,"Position"],truncation=c(0,20),weighting="1")
mca1 <- mca(y=log(Salmon[,"Abundance"]+1),x=Salmon[,"Substrate"],memobj=map)
mca1
plot(mca1)
summary(mca1)   # Works only with tested mca.
#
## Parametric test.
mcatest1 <- test.mca(mca1)
mcatest1
summary(mcatest1)
plot(mcatest1)
plot(Abundance~Position,data=Salmon, type="l")
lines(y=exp(fitted(mcatest1)-1),x=Salmon[,"Position"])
plot(x=Salmon[,"Abundance"],y=exp(fitted(mcatest1)-1),asp=1,xlim=c(0,14),ylim=c(0,14))
abline(0,1)
residuals(mcatest1)
#
## Permutation test. Warning: takes time to calculate.
mcapermute1 <- permute.mca(mca1)
mcapermute1
summary(mcapermute1)
plot(Abundance~Position,data=Salmon, type="l")
lines(y=exp(fitted(mcapermute1)-1),x=Salmon[,"Position"])
plot(x=Salmon[,"Abundance"],y=exp(fitted(mcapermute1)-1),asp=1,xlim=c(0,14),ylim=c(0,14),
     xlab="Observed abundance",ylab="Fitted abundance")
abline(0,1)
residuals(mcapermute1)
#
\dontrun{
### Exemple 2: Oribatid mites in Lake Cromwell, Quebec, Canada.
#
## Requires package vegan
library(vegan)
data(mite.xy)
map <- meigmap(x=as.matrix(mite.xy),truncation=c(0,NA),weighting="f1")
data(mite)
data(mite.env)
mca2 <- mca(y=log(mite[,"LRUG"]+1),x=mite.env[,"WatrCont"],memobj=map)
mca2
mcapermute2 <- permute.mca(mca2)
summary(mcapermute2)
#
layout(matrix(1:2,1,2))
obs <- 0.5+log(mite[,"LRUG"]+1)/2
plot(y~x,data=mite.xy, asp=1, pch=21, cex=obs, bg="black",main="Observed")
fit <- 0.5+fitted(mcapermute2)/2
plot(y~x,data=mite.xy, asp=1, pch=21, cex=fit, bg="black",main="Fitted")
#
layout(1)
plot(x=log(mite[,"LRUG"]+1),y=fitted(mcapermute2),asp=1,xlim=c(0,7),ylim=c(0,7),
     xlab="Observed abundance (log(x+1))",ylab="Fitted abundance (log(x+1))")
abline(0,1)
#
}
}
