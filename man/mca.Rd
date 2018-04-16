\encoding{utf8}
\name{MCA}
\alias{MCA}
\alias{test.cdp}
\alias{permute.cdp}
\alias{parPermute.cdp}
\title{ Multiple-descriptors, Multiscale Codependence Analysis }
\description{ Functions to perform Multiscale Codependence Analysis
  (MCA) }
\usage{
MCA(Y, X, emobj)
test.cdp(object, alpha = 0.05, max.step, response.tests = TRUE)
permute.cdp(object, permute, alpha = 0.05, max.step,
            response.tests = TRUE)
parPermute.cdp(object, permute, alpha = 0.05, max.step,
               response.tests = TRUE, nnode, seeds, verbose = TRUE,
               ...)
}
\arguments{
  \item{Y}{ a numeric matrix or vector containing the response
    variable(s). }
  \item{X}{ a numeric matrix or vector containing the explanatory
    variable(s). }
  \item{emobj}{ a \code{\link{eigenmap-class}} object. }
  \item{object}{ a \code{\link{cdp-class}} object. }
  \item{alpha}{ type I (\eqn{\alpha}{alpha}) error threshold used by the
    testing procedure. }
  \item{max.step}{ maximum number of steps to perform when testing for
    statistical significance. }
  \item{response.tests}{ a boolean specifying whether to test individual
    response variables. }
  \item{permute}{ The number of random permutations used for
    testing. When omitted, the number of permutations is calculated
    using function \code{\link{minpermute}}.}
  \item{nnode}{ The number of parallel computation nodes. }
  \item{seeds}{ Random number generator seeds for parallel the
    computation nodes. }
  \item{verbose}{ Whether to return user notifications. }
  \item{...}{ Parameters to be passed to \code{parallel::makeCluster()}}
}
\details{
  Multiscale Codependence Analysis (MCA) allows to calculate
  correlation-like (i.e.codependence) coefficients between two variables
  with respect to structuring variables (Moran's eigenvector maps). The
  purpose of this function is limited to parameter fitting. Test
  procedures are handled through \code{test.cdp} (parametric testing)
  or \code{permute.cdp} (permutation testing). Additionaly, methods are
  provided for printing, displaying the testing summary, plotting
  results, calculating fitted and residuals values, and making
  predictions. It is noteworthy that the test procedure used by
  \code{MCA} deviates from the standard R workflow since intermediate
  testing functions (\code{test.cdp} and \code{permute.cdp}) need
  first to be called before any testing be performed. For \code{MCA},
  testing functionalities had been moved away from summary.cdp because
  testing is computationally intensive. Function \code{parPermute.cdp}
  allows the user to spread the number of permutation on many
  computation nodes. It relies on package \code{parallel}. Omitting
  parameter \code{nnode} lets function \code{parallel::detectCores()}
  specify the number of node. Similarly, omitting parameter \code{seeds}
  lets the draw seeds uniformly between
  \code{±.Machine$integer.max}. If needed, one may pass initialization
  parameters to \code{parallel::makeCluster()}.
}
\value{
  A \code{\link{cdp-class}} object.
}
\references{
  Guénard, G., Legendre, P., Boisclair, D., and Bilodeau, M. 2010.
  Multiscale codependence analysis: an integrated approach to analyse
  relationships across scales. Ecology 91: 2952-2964
  
  Guénard, G. Legendre, P. 2018. Bringing multivariate support to
  multiscale codependence analysis: Assessing the drivers of community
  structure across spatial scales. Meth. Ecol. Evol. 9: 292-304
}
\author{ Guillaume Guénard, Département des sciences biologiques,
  Université de Montréal, Montréal, Québec, Canada. }
\seealso{
  \code{\link{eigenmap}}
}
\examples{
#
###### Begin {Salmon exemple}
#
data(Salmon)
#
## Converting the data from data frames to to matrices:
Abundance <- log1p(as.matrix(Salmon[,"Abundance",drop=FALSE]))
Environ <- as.matrix(Salmon[,3L:5])
#
## Creating a spatial eigenvector map:
map1 <- eigenmap(x=Salmon[,"Position"],weighting=Wf.binary,boundaries=c(0,20))
#
## Case of a single descriptor:
mca1 <- MCA(Y=Abundance,X=Environ[,"Substrate",drop=FALSE],emobj=map1)
mca1
mca1_partest <- test.cdp(mca1)
mca1_partest
summary(mca1_partest)
par(mar = c(6,4,2,4))
plot(mca1_partest, las = 3)
mca1_pertest <- permute.cdp(mca1)
\dontrun{
## or:
mca1_pertest <- parPermute.cdp(mca1,permute=999999)
}
mca1_pertest
summary(mca1_pertest)
plot(mca1_pertest, las = 3)
mca1_pertest$UpYXcb$C # Array containing the codependence coefficients
#
## With all descriptors at once:
mca2 <- MCA(Y=log1p(as.matrix(Salmon[,"Abundance",drop=FALSE])),
            X=as.matrix(Salmon[,3L:5]),emobj=map1)
mca2
mca2_partest <- test.cdp(mca2)
mca2_partest
summary(mca2_partest)
par(mar = c(6,4,2,4))
plot(mca2_partest, las = 3)
mca2_pertest <- permute.cdp(mca2)
\dontrun{
or:
mca2_pertest <- parPermute.cdp(mca2,permute=999999)
}
mca2_pertest
summary(mca2_pertest)
plot(mca2_pertest, las = 3)
mca2_pertest$UpYXcb$C # Array containing the codependence coefficients
mca2_pertest$UpYXcb$C[,1L,] # now turned into a matrix.
#
###### End {Salmon exemple}
#
###### Begin {Doubs exemple}
#
data(Doubs)
#
## Creating a spatial eigenvector map:
map2 <- eigenmap(x=Doubs.geo[,"DFS"])
#
mca3 <- MCA(Y=log1p(Doubs.fish),X=Doubs.env,emobj=map2)
mca3
mca3_pertest <- permute.cdp(mca3)
\dontrun{
## or:
mca3_pertest <- parPermute.cdp(mca3,permute=999999)
}
mca3_pertest
summary(mca3_pertest)
par(mar = c(6,4,2,4))
plot(mca3_pertest, las = 2)
mca3_pertest$UpYXcb$C # Array containing the codependence coefficients
#
## Display the results along the transect
spmeans <- colMeans(log1p(Doubs.fish))
pca1 <- svd(log1p(Doubs.fish) - rep(spmeans,each=nrow(Doubs.fish)))
par(mar = c(5,5,2,5)+0.1)
plot(y = pca1$u[,1L], x = Doubs.geo[,"DFS"], pch = 21L, bg = "red",
     ylab = "PCA1 loadings", xlab = "Distance from river source (km)")
#
x <- seq(0,450,1)
newdists <- matrix(NA, length(x), nrow(Doubs.geo))
for(i in 1L:nrow(newdists))
  newdists[i,] <- abs(Doubs.geo[,"DFS"] - x[i])
#
## Calculating predictions for arbitrary sites under the same set of
## environmental conditions that the codependence model was built with.
prd1 <- predict(mca3_pertest,
                newdata=list(target = eigenmap.score(map2, newdists)))
#
## Projection of the predicted species abundance on pca1:
Uprd1 <- (prd1 - rep(spmeans, each = nrow(prd1))) \%*\% pca1$v \%*\% diag(pca1$d^-1)
lines(y = Uprd1[,1L], x = x, col=2, lty = 1)
#
## Projection of the predicted species abundance on pca2:
plot(y = pca1$u[,2L], x = Doubs.geo[,"DFS"], pch = 21L, bg = "red",
     ylab = "PCA2 loadings", xlab = "Distance from river source (km)")
lines(y = Uprd1[,2L], x = x, col=2, lty = 1)

#
## Displaying only the observed and predicted abundance for Brown Trout.
par(new=TRUE)
plot(y = log1p(Doubs.fish[,"TRU"]),Doubs.geo[,"DFS"],pch=21L,bg="green",
     ylab="",xlab="",new=FALSE,axes=FALSE)
axis(4)
lines(y = prd1[,"TRU"], x = x, col=3)
mtext(side=4, "log(Abundance+1)", line = 2.5)
#
###### End {Doubs exemple}
#
###### Begin {Oribatid exemple}
#
data(Mite)
#
map3 <- eigenmap(x = mite.geo)
# Organize the environmental variables
mca4 <- MCA(Y = log1p(mite.species), X = mite.env, emobj = map3)
mca4_partest <- test.cdp(mca4, response.tests = FALSE)
summary(mca4_partest)
plot(mca4_partest, las = 2, lwd = 2)
plot(mca4_partest, col = rainbow(1200)[1L:1000], las = 3, lwd = 4,
     main = "Codependence diagram", col.signif = "white")
#
rng <- list(x = seq(min(mite.geo[,"x"]) - 0.1, max(mite.geo[,"x"]) + 0.1, 0.05),
            y = seq(min(mite.geo[,"y"]) - 0.1, max(mite.geo[,"y"]) + 0.1, 0.05))
grid <- cbind(x = rep(rng[["x"]], length(rng[["y"]])),
              y = rep(rng[["y"]], each = length(rng[["x"]])))
newdists <- matrix(NA, nrow(grid), nrow(mite.geo))
for(i in 1L:nrow(grid)) {
  newdists[i,] <- ((mite.geo[,"x"] - grid[i,"x"])^2 +
                    (mite.geo[,"y"] - grid[i,"y"])^2)^0.5
}
#
spmeans <- colMeans(mite.species)
pca2 <- svd(log1p(mite.species) - rep(spmeans, each = nrow(mite.species)))
#
prd2 <- predict(mca4_partest,
          newdata = list(target = eigenmap.score(map3, newdists)))
Uprd2 <- (prd2 - rep(spmeans, each = nrow(prd2))) \%*\% pca2$v \%*\% diag(pca2$d^-1)
#
### Printing the response variable
prmat <- Uprd2[,1L]
dim(prmat) <- c(length(rng$x),length(rng$y))
zlim <- c(min(min(prmat),min(pca2$u[,1L])),max(max(prmat),max(pca2$u[,1L])))
image(z = prmat, x = rng$x, y = rng$y, asp = 1, zlim = zlim,
      col = rainbow(1200L)[1L:1000], ylab = "y", xlab = "x")
points(x = mite.geo[,"x"], y = mite.geo[,"y"], pch = 21,
  bg = rainbow(1200L)[round(1+(999*(pca2$u[,1L]-zlim[1L])/(zlim[2L]-zlim[1L])),0)])
#
###### End {Oribatid exemple}
#
}
