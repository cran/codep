% --- Source file: ./man/mevmap.Rd ---
% File src/library/codep/man/mevmap.Rd
% Part of the R package, http://www.R-project.org
% Copyright 2009 Guenard, Guillaume and Legendre, Pierre
% Distributed under GPL 2 or later
\name{meigmap}
\alias{meigmap}
\alias{mem}
\alias{print.mem}
\alias{plot.mem}
\title{Moran's eigenvector maps}
\description{
  Functions and methods to handle Moran's eigenvector maps of a set of
  locations in a space with an arbitrary number of dimension.
}
\usage{
meigmap(x, opt.coordinates=NA, truncation=c(0,NA),weighting=c("1","f1","f2","f3"), wpar=1, select=1e-9)
mem(coordinates, lambda, U)
\method{print}{mem}(x, ...)
\method{plot}{mem}(x, ...)
}
\arguments{
  \item{x}{ A set of coordinates defined in one (numeric vector) or many
    (a coordinate x dimension matrix) dimensions or, alternatively, a
    distance matrix provided by \code{\link{dist}}. Coordinates are
    treated as geographic and distances between them are taken to be
    Euclidean. For \code{print} and \code{plot} methods: an object of
    \code{\link{class}} \sQuote{mem}. }
  \item{opt.coordinates}{ Coordinates to be used when a distance matrix
    is provided as x. Used for plotting purposes. }
  \item{truncation}{ Threshold values (minimum and maximum) used to
    obtain the connectivity matrix. The minimum value is the value above
    which locations are to be considered as neighbour (in most case 0)
    whereas the maximum value is the value above which locations are to
    be considered as distant. If any trunctation value is \code{NA}, the
    minimum value is taken to be 0 and the maximum value is taken as the
    minimum value that allow every locations to be either directly or
    indirectly connected with one another. }
  \item{weighting}{ The function used to obtain the weighting matrix:
    
    1 refers to the identity function and implies that the spatial
    weighting matrix is equal to the connnectivity matrix,
    
    f1 is \eqn{d_{i,j} = 1 - (d_{i,j} / d_max)}{d i,j = 1 - (d i,j /
    dmax)} where \eqn{d_max}{dmax} is the distance between the two most
    distant locations in theset,
  
    f2 is \eqn{d_{i,j} = 1 - (d_{i,j} / d_max)^wpar}{d i,j = 1 - (d i,j
    / dmax)^wpar}, and
    
    f3 is \eqn{d_{i,j} = 1 / d_{i,j}^wpar}{d i,j = 1 / (d i,j)^wpar}.
  }
  \item{wpar}{ Parameter of the weighting function f2 or f3. Ignored
    when weighting is 1 or f1. }
  \item{select}{ The smallest absolute eigenvalue for eigenfunctions to
    be considered as a suitable variables. }
  \item{coordinates}{ A set of n coordinates for the mem object. }
  \item{lambda}{ The p eigenvalues of the mem object. }
  \item{U}{ A n x p matrix of structuring variables. }
  \item{...}{ Further parameters to be passed to other functions or
    methods (currently ignored). }
}
\details{
  Moran's eigenvector maps are sets of orthonormal functions
  (i.e. eigenfunction) obtained from the locations of observations
  within a structuring framework, e.g., space, time, along graph. It is
  obtained by the eigenvalue decomposition of a spatial weighting
  matrix; a centred square matrix itself obtained from a distance
  matrix. The spatial weighting matrix is the haddamar product of a
  connectivity matrix \eqn{B}{B} and a weighting matrix \eqn{A}{A}. The
  function described herein handles user-defined truncations parameters
  to calculate \eqn{B}{B} and provides a default approach to estimate
  these parameters should they be missing. It also offers four different
  ways to compute \eqn{A}{A} through parameters \code{weighting} and
  \code{wpar}.
  
  Function \code{mem} defines a mem object customly, by providing
  coordinates, eigenvalues, and structuring variables as column
  vectors. Basin checks for size compatibility are made. It falls on
  user's responsability to further ensure the adequacy of the object.
  
  The \code{print} method provides the number of the number of
  orthonormal variables (i.e. basis functions), the number of
  observations these functions are spanning, and their associated
  eigenvalues.
  
  The \code{plot} method provides a plot of the eigenvalues and offers
  the possibility to plot the values of variables for 1- or
  2-dimensional sets of coordinates. To use this later functionnality, a
  coordinate parameter must be provided with the actual plotting
  coordinates. \code{plot.mem} opens the default graphical device
  driver, i.e., \code{X11}, \code{windows}, or \code{quartz} and
  recurses through variable with a left mouse click on the graphical
  window. A right mouse click interrupts recursing on \code{X11} and
  \code{windows}.
  Recursing is stopped by hitting \emph{Esc} on the \code{quartz}
  graphical device driver (Mac OS X users).
  
  See Dray et al. (2006) for further details about this method.
}
\value{
  An object of \code{\link{class}} \sQuote{mem} is a list containing:
  \item{coordinates}{A matrix of coordinates associated with }
  \item{lambda}{ A vector of the eigenvalues obtain from the computation
    of the the Moran's eigenvector map. }
  \item{U}{ A matrix of the eigenvectors defining the Moran's
    eigenvector map. }
}
\references{
Dray, S.; Legendre, P. and Peres-Neto, P. 2006. Spatial modelling: a
comprehensive framework for principal coordinate analysis of neighbor
matrices (PCNM). Ecol. Modelling 196: 483-493
}

\author{ Guillaume Guenard, Laboratoire evolution et diversite
  biologique, CNRS / Universite Paul-Sabatier, Toulouse, France.
}

\seealso{
  \code{\link{mca}}
}
\examples{
#
### Exemple 1: A linear transect.
#
data(Salmon)
map <- meigmap(x=Salmon[,"Position"],truncation=c(0,20),weighting="1")
map
plot(map)
#
### Exemple 2: A unevenly sampled surface.
#
## Require package vegan
library(vegan)
data(mite.xy)
map <- meigmap(x=as.matrix(mite.xy),truncation=c(0,NA),weighting="f1")
map
plot(map)
#
### Exemple 3: Building a mem object for other basis functions (fourier transforms).
#
U <- matrix(NA,100,0)
for (i in 1:49) {
  for (j in c(1+0i, 0+1i)) {
    spc <- complex(100)
    spc[i+1] <- j
    spc[100-i+1] <- Conj(j)
    U <- cbind(U,Re(fft(spc,inverse=TRUE))/sqrt(200))
  }
}
spc <- complex(100)
spc[51] <- 1+0i
U <- cbind(U,Re(fft(spc,inverse=TRUE))/10)
dmat <- matrix(0,100,100)
for (i in 1:99) {
  dmat[i,i+1] <- 1
  dmat[i+1,i] <- 1
}
dmat[1,100] <- dmat[100,1] <- 1
dmat <- (diag(100)-(rep(0.1,100)\%*\%t(rep(0.1,100))))\%*\%dmat\%*\%(diag(100)-(rep(0.1,100)\%*\%t(rep(0.1,100))))
lambda <- diag(t(U) \%*\% dmat \%*\% U)
map <- mem(coordinates=1:100, lambda=lambda, U=U)
map
plot(map)
#
}
