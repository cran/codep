\encoding{utf8}
\name{MCA}
\alias{MCA}
\alias{test.mca}
\alias{permute.mca}
\title{ Multiscale Codependence Analysis }
\description{ Functions to perform Multiscale Codependence Analysis
  (MCA) }
\usage{
MCA(y, x, emobj)
test.mca(mcaobj, alpha = 0.05, max.step = Inf)
permute.mca(mcaobj, permute = NA, alpha = 0.05, max.step = Inf)
}
\arguments{
  \item{y}{ A numeric vector containing the response variable. }
  \item{x}{ A numeric vector containing the explanatory variable for
    functions \code{test.mca} or \code{permute.mca}. }
  \item{emobj}{ A \code{\link{eigenmap-class}} object. }
  \item{mcaobj}{ A \code{\link{mca-class}} object. }
  \item{alpha}{ The type I (\eqn{\alpha}{alpha}) error threshold used by
    the testing procedure. }
  \item{max.step}{ The maximum number of steps to perform when testing
    for statistical significance. }
  \item{permute}{ The number of random permutations used for
    testing. Obtained through \code{\link{minpermute}} when let to
    \code{NA}. }
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
  values, and making predictions. It is noteworthy that the test
  procedure used by \code{\link{MCA}} deviates from the standard R
  workflow since intermediate testing functions (\code{\link{test.mca}}
  and \code{\link{permute.mca}}) need first to be called before any
  testing be performed. For \code{\link{MCA}} testing functionalities
  had been moved away from summary.mca because of permutation testing,
  which is the most relevant option in practice, takes long time to
  run.}
\value{
  A \code{\link{mca-class}} object.
}
\references{
  Guénard, G., Legendre, P., Boisclair, D., and Bilodeau, M. 2010.
  Multiscale codependence analysis: an integrated approach to analyse
  relationships across scales. Ecology 91: 2952-2964
}
\author{ Guillaume Guénard, Département des sciences biologiques,
  Université de Montréal, Montréal, Québec, Canada. }
\seealso{
  \code{\link{eigenmap}}
}
\examples{
#
### Example 1: Atlantic salmon in the St-Marguerite river, Quebec, Canada.
#
data(Salmon)
map <- eigenmap(x=Salmon[,"Position"],truncate=c(0,20),wf="binary")
mca1 <- MCA(y=log(Salmon[,"Abundance"]+1),x=Salmon[,"Substrate"],emobj=map)
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
### Example 2: Oribatid mites in Lake Cromwell, Quebec, Canada.
#
## Requires package vegan
library(vegan)
data(mite.xy)
map <- eigenmap(x=as.matrix(mite.xy),truncate=c(0,NA),wf="Drayf1")
data(mite)
data(mite.env)
mca2 <- MCA(y=log(mite[,"LRUG"]+1),x=mite.env[,"WatrCont"],emobj=map)
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
