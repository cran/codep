\encoding{utf8}
\name{eigenmap}
\alias{eigenmap}
\title{Spatial eigenvector maps}
\description{
  Function to spatial eigenvector maps of a set of locations in a space
  with an arbitrary number of dimension.
}
\usage{
eigenmap(x,opt.coord=NA,truncate=c(0,NA),wf,wpar,select=.Machine$double.eps^0.5)
}
\arguments{
  \item{x}{ A set of coordinates defined in one (numeric vector) or many
    (a coordinate x dimension matrix) dimensions or, alternatively, a
    distance matrix provided by \code{\link{dist}}. Coordinates are
    treated as geographic and distances between them are taken to be
    Euclidean. }
  \item{opt.coord}{ Coordinates to be used when a distance matrix is
    provided as x. Used for plotting purposes. }
  \item{truncate}{ Threshold values (minimum and maximum) used to obtain
    the connectivity matrix. Pairs of location whose distance to one
    another are between these values are considered as neighbours
    (\eqn{b_{i,j}=1}{b_{i,j}=1}) whereas values located below the
    mininum and above the maximum are considered as equivallent or
    distant, respectively (\eqn{b_{i,j}=0}{b_{i,j}=0} in both
    cases). Defaults are \eqn{0}{0} for the minimum value and \code{NA}
    for the maximum. Values \code{NA} indicates the function to take the
    minimum value that allow every locations to form a single cluster
    following single linkage clustering as a maximum value (obtained
    from \code{\link{hclust}}}. Ignored when \code{weighting="sqrd"}.
  \item{wf}{ Choice of a function to obtain the edge weighting matrix:
    
    binary (default value) refers to the identity function (The spatial
    weighting matrix is simply the connectivity matrix),

    \code{dbMEM} is \eqn{a_{i,j} = 1 - (d_{i,j} /
      (4*truncate_2))^2}{a_{i,j} = 1 - (d_{i,j} / (4*truncate_2))^2}
    
    \code{Drayf1} is \eqn{a_{i,j} = 1 - (d_{i,j} / d_{max})}{a_{i,j} = 1 -
      (d_{i,j} / d_max)} where \eqn{d_max}{dmax} is the distance between
      the two most distant locations in the set,
   
    \code{Drayf2} is \eqn{a_{i,j} = 1 - (d_{i,j} / d_{max})^{wpar}}{a_{i,j}
      = 1 - (d_{i,j} / d_max)^wpar},
    
    \code{Drayf3} is \eqn{a_{i,j} = 1 / d_{i,j}^{wpar}}{a_{i,j} = 1 /
     d_{i,j}^wpar}, and

    \code{sqrd} consists in taking \eqn{w_{i,j} = -0.5*d_{i,j}}{w_{i,j}
      = -0.5*d_{i,j}} and does not involve any trunction.

    Functions \code{f1}, \code{f2}, and \code{f3} were proposed by Dray
    et al. (2006) and function \code{dbMEM} was proposed by Legendre and
    Legendre (2013).

    The \code{sqrd} weighting approach is equivallent to submitting the
    elementwise square-root of the distance matrix to a principal
    coordinate analysis. That option is not much documented in the
    ecological litterature, but is actually equivallent, for evenly
    spaced transect or surfaces (square or rectangle), to using the
    basis functions of type II discrete cosine basis transforms. }
  \item{wpar}{ Exponent in weighting function \code{Drayf2} and
    \code{Drayf3}. Ignored when weighting is \code{binary},
    \code{dbMEM}, \code{Drayf1}, or \code{sqrd}. }
  \item{select}{ The smallest absolute eigenvalue for eigenfunctions to
    be considered as a suitable predictive variables. Default value
    depends on one's particular machine and is set to the square-root of
    .Machine$double.eps }
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
}
\value{
  A \code{\link{eigenmap-class}} object.
}
\references{
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
map <- eigenmap(x=Salmon[,"Position"],truncate=c(0,20),wf="sqrd")
map
plot(map)   # Same as basis functions in type II Discrete Cosine Transforms.
#
map <- eigenmap(x=Salmon[,"Position"],truncate=c(0,20),wf="Drayf1")
map
plot(map)   # Seem like some undifined-type Discrete Sine Transform basis functions
#
map <- eigenmap(x=Salmon[,"Position"],truncate=c(0,20),wf="Drayf2",wpar=2)
map
plot(map)
#
map <- eigenmap(x=Salmon[,"Position"],truncate=c(0,20),wf="Drayf3",wpar=2.4)
map
plot(map)
#
map <- eigenmap(x=Salmon[,"Position"],truncate=c(0,20),wf="dbMEM")
map
plot(map)   # Similar to Drayf1.
#
### Example 2: A unevenly sampled surface.
#
## Require package vegan
\dontrun{
library(vegan)
data(mite.xy)
map <- eigenmap(x=as.matrix(mite.xy),truncate=c(0,NA),wf="sqrd")
map
plot(map)
}
}
