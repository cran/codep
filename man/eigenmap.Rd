\encoding{utf8}
\name{eigenmap}
\alias{eigenmap}
\alias{Wf.sqrd}
\alias{Wf.RBF}
\alias{Wf.binary}
\alias{Wf.PCNM}
\alias{Wf.Drayf1}
\alias{Wf.Drayf2}
\alias{Wf.Drayf3}
\alias{eigenmap.score}
\title{Spatial eigenvector maps}
\description{
  Function to calculate spatial eigenvector maps of a set of locations
  in a space with an arbitrary number of dimension.
}
\usage{
eigenmap(x,opt.coord=NA,weighting=Wf.sqrd,boundaries,wpar,select=.Machine$double.eps^0.5)
Wf.sqrd(D)
Wf.RBF(D,wpar=1)
Wf.binary(D,boundaries)
Wf.PCNM(D,boundaries)
Wf.Drayf1(D,boundaries)
Wf.Drayf2(D,boundaries,wpar=1)
Wf.Drayf3(D,boundaries,wpar=1)
eigenmap.score(object,target)
}
\arguments{
  \item{x}{ A set of coordinates defined in one (numeric vector) or many
    (a coordinate x dimension matrix) dimensions or, alternatively, a
    distance matrix provided by \code{\link{dist}}. Coordinates are
    treated as cartesian coordinates and the distances between them are
    assumed to be Euclidean. }
  \item{opt.coord}{ Coordinates to be used when a distance matrix is
    provided as x. Used for plotting purposes. }
  \item{weighting}{ The function to obtain the edge weighting
    matrix. That function must have the raw distances as a first
    parameter, optionally a second parameter named \code{boundaries}
    giving the boundaries of the within which locations are regarded as
    neighbour and a third parameter named \code{wpar} containing any
    other weighting function parameter.

    \code{Wf.sqrd} consists in taking \eqn{w_{i,j} =
      -0.5*d_{i,j}}{w_{i,j} = -0.5*d_{i,j}} and does not involve any
    trunction.

    \code{Wf.sqrd} consists in taking \eqn{w_{i,j} =
      exp(-wpar*d_{i,j}^2)}{w_{i,j} = exp(-wpar*d_{i,j}^2)} and does not
    involve any trunction.

    \code{Wf.binary} (default value) the spatial weighting matrix is
    simply the connectivity matrix,

    \code{Wf.PCNM} is \eqn{a_{i,j} = 1 - (d_{i,j} /
      (4*boundaries_2))^2}{a_{i,j} = 1 - (d_{i,j} / (4*boundaries_2))^2}

    \code{Wf.Drayf1} is \eqn{a_{i,j} = 1 - (d_{i,j} / d_{max})}{a_{i,j}
      = 1 - (d_{i,j} / d_max)} where \eqn{d_max}{dmax} is the distance
      between the two most distant locations in the set,

    \code{Wf.Drayf2} is \eqn{a_{i,j} = 1 - (d_{i,j} /
      d_{max})^{wpar}}{a_{i,j} = 1 - (d_{i,j} / d_max)^wpar},

    \code{Wf.Drayf3} is \eqn{a_{i,j} = 1 / d_{i,j}^{wpar}}{a_{i,j} = 1 /
     d_{i,j}^wpar}, and

    Functions \code{Wf.Drayf1}, \code{Wf.Drayf2}, and \code{Wf.Drayf3}
    were proposed by Dray et al. (2006) and function \code{PCNM} was
    proposed by Legendre and Legendre (2012).
    The \code{Wf.sqrd} weighting approach is equivallent to submitting
    the elementwise square-root of the distance matrix to a principal
    coordinate analysis. That option is not much documented in the
    ecological litterature, but is actually equivallent, for evenly
    spaced transect or surfaces (square or rectangle), to using the
    basis functions of type II discrete cosine basis transforms. }
   \item{boundaries}{ (optional) Threshold values (minimum and maximum)
    used to obtain the connectivity matrix. Pairs of location whose
    distance to one another are between these values are considered as
    neighbours (\eqn{b_{i,j}=1}{b_{i,j}=1}) whereas values located below
    the mininum and above the maximum are considered as equivallent or
    distant, respectively (\eqn{b_{i,j}=0}{b_{i,j}=0} in both
    cases). Defaults are \eqn{0}{0} for the minimum value and \code{NA}
    for the maximum. Values \code{NA} indicates the function to take the
    minimum value that allow every locations to form a single cluster
    following single linkage clustering as a maximum value (obtained
    from \code{\link{hclust}}}. Ignored when \code{weighting="Wf.sqrd"}
    or \code{weighting="Wf.RBF"}.
  \item{wpar}{ Weighting function parameters. }
  \item{select}{ The smallest absolute eigenvalue for eigenfunctions to
    be considered as a suitable predictive variables. Default value
    depends on one's particular computer and is set to the square-root
    of \code{.Machine$double.eps} }
  \item{D}{ A distance matrix. }
  \item{object}{ A \code{\link{eigenmap-class}} object. }
  \item{target}{ A set of distances between the sampling locations
    (passed to \code{\link{eigenmap}} using \code{x}) and the target
    locations where spatially-explicit predictions are to be made. }
}
\details{
  Spatial eigenvector maps are sets of eigenfunctions obtained from the
    locations of the observations in a structuring framework, e.g.,
    space, time, or in a graph. It is obtained by eigenvalue
    decomposition of a  spatial weighting matrix, computed as described
    in Dray et al. (2006) and Legendre & Legendre (2012, Section
    14.2). That square matrix is Gower-centred before
    eigen-decomposition. The spatial weighting matrix is the Hadamard
    product of a connectivity matrix \eqn{\mathbf{B}}{B} and an edge
    weighting matrix \eqn{\mathbf{A}}{A}. The function described herein
    handles user-chosen truncation parameters to calculate
    \eqn{\mathbf{B}}{B} and provides a default approach to estimate
    these parameters should they be missing. It also offers four
    different ways of computing \eqn{\mathbf{A}}{A} through parameters
    \code{weighting} and \code{wpar}.

    In is noteworthy that in the present implementation, matrix
    \eqn{\mathbf{B}}{B} is not obtained using a minimum spanning tree as
    suggested by Dray et al. (2006) but using a simpler approach whereby
    every distances within a user-defined trunction interval are taken
    as neighbour.

    Functions \code{\link{Wf.sqrd}}, \code{\link{Wf.RBF}},
    \code{\link{Wf.binary}}, \code{\link{Wf.PCNM}},
    \code{\link{Wf.Drayf1}}, \code{\link{Wf.Drayf2}}, and
    \code{\link{Wf.Drayf3}} are not intested to be called as is but
    through \code{\link{eigenmap}} (and within
    \code{\link{eigenmap.score}}). Other, user-defined, function can be
    used by \code{\link{eigenmap}} and should be visible to if one wants
    to call \code{\link{eigenmap.score}} to obtain predictors.

    For \code{\link{eigenmap.score}}, the distances between sampling
    locations and the targets locations must be of the same type as
    those that had been passed to \code{\link{eigenmap}}. If cartesian
    coordinates were passed to \code{x}, the distances to target must be
    Euclidean.
}
\value{
  \code{\link{eigenmap}} returns a \code{\link{eigenmap-class}} object
    and \code{\link{eigenmap.score}} returns a the scores on for each
    target locations 
}
\references{
  Borcard, D. and Legendre, P. 2002. All-scale spatial analysis of
  ecological data by means of principal coordinates of neighbour
  matrices. Ecol. Model. 153: 51-68
  
  Dray, S.; Legendre, P. and Peres-Neto, P. 2006. Spatial modelling: a
  comprehensive framework for principal coordinate analysis of neighbor
  matrices (PCNM). Ecol. Modelling 196: 483-493

  Legendre, P. and Legendre, L. 2012. Numerical Ecology, 3rd English
  edition. Elsevier Science B.V., Amsterdam, The Neatherlands.
}
\author{ Guillaume Guénard, Departement des sciences biologiques,
  Universite de Montréal, Montréal, Quebec, Canada. }

\seealso{
  \code{\link{MCA}}
  \code{\link{eigenmap-class}}
}
\examples{
#
### Example 1: A linear transect.
#
data(Salmon)
#
## No boundaries provided for a function that requires them: a warning is issued
map <- eigenmap(x=Salmon[,"Position"],weighting=Wf.binary)
map # plot(map)
#
## Boundaries are provided: the function is happy
map <- eigenmap(x=Salmon[,"Position"],weighting=Wf.binary,boundaries=c(0,20))
map # plot(map)
#
map <- eigenmap(x=Salmon[,"Position"],weighting=Wf.Drayf1,boundaries=c(0,20))
map # plot(map)
#
map <- eigenmap(x=Salmon[,"Position"],weighting=Wf.Drayf2,boundaries=c(0,20))
map # plot(map)
#
map <- eigenmap(x=Salmon[,"Position"],weighting=Wf.Drayf3,boundaries=c(0,20),wpar=2)
map # plot(map)
#
map <- eigenmap(x=Salmon[,"Position"],weighting=Wf.PCNM,boundaries=c(0,20))
map # plot(map)
#
map <- eigenmap(x=Salmon[,"Position"],weighting=Wf.sqrd)
map # plot(map)
#
map <- eigenmap(x=Salmon[,"Position"],weighting=Wf.RBF,wpar=0.001)
map # plot(map)
#
### Example 2: Using predictior scores
#
smpl <- c(4,7,10,14,34,56,61,64)  # A sample to discard
map <- eigenmap(x=Salmon[-smpl,"Position"],weighting=Wf.sqrd)
scr <- eigenmap.score(object=map,target=as.matrix(dist(Salmon[,"Position"]))[,-smpl])
all(round(scr[-smpl,] - map$U, 10) == 0) # Scores of sampling points are the eigenvectors
scr[smpl,]
#
wh <- 5L   # You can try with other vectors.
plot(map$U[,wh]~Salmon[-smpl,"Position"], ylab = expression(U[5]),
     xlab = "Position along transect")
points(y=scr[smpl,wh],x=Salmon[smpl,"Position"],pch=21,bg="black")
#
map <- eigenmap(x=Salmon[-smpl,"Position"],weighting=Wf.binary,boundaries=c(0,20))
scr <- eigenmap.score(object=map,target=as.matrix(dist(Salmon[,"Position"]))[smpl,-smpl])
#
wh <- 1L   # You can try with other vectors.
plot(map$U[,wh]~Salmon[-smpl,"Position"], ylab = expression(U[1]),
     xlab = "Position along transect (m)")
points(y=scr[,wh],x=Salmon[smpl,"Position"],pch=21,bg="black")
#
map <- eigenmap(x=Salmon[-smpl,"Position"],weighting=Wf.PCNM,boundaries=c(0,100))
scr <- eigenmap.score(object=map,target=as.matrix(dist(Salmon[,"Position"]))[smpl,-smpl])
#
wh <- 1L   # You can try with other vectors.
plot(map$U[,wh]~Salmon[-smpl,"Position"], ylab = expression(U[1]),
     xlab = "Position along transect (m)")
points(y=scr[,wh],x=Salmon[smpl,"Position"],pch=21,bg="black")
#
### Example 3: A unevenly sampled surface.
#
data(Mite)
map <- eigenmap(x=as.matrix(mite.geo),weighting=Wf.sqrd)
map # plot(map)
#
map <- eigenmap(x=as.matrix(mite.geo),weighting=Wf.RBF)
map # plot(map)
#
}
