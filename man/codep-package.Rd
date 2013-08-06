\encoding{utf8}
\name{codep-package}
\alias{codep-package}
\docType{package}
\title{ Functions and methods to perform Multiscale Codependence
  Analysis (MCA) }
\description{
  Multiscale Codependence Analysis (MCA) consists in assessing the
  coherence of pairs of variables is space (or time) using the product
  of their correlation coefficients with series of spatial (or temporal)
  eigenfunctions. That product, which is positive or negative when
  variables show similar or opposing trends, respectively, are called
  codependence coefficients. These eigenfunctions are obtained in three
  steps: 1) a distance matrice calculated from the locations of samples
  in space (or the organisation of the sampling schedule). 2) from that
  distance matrix, a matrix of spatial weights is obtained; the same
  matrix as to calculate Moran's autocorrelation index, hence the name,
  and 3) the spatial weight matrix is eigenvalue-decomposed after
  centering each rows and columns of the spatial weight matrix.

  The statistical significance of the codependence coefficients is
  tested using parametric or permutational testing of a \eqn{\tau}{tau}
  statistic. The \eqn{\tau}{tau} statistic is the product of the two
  Student's \eqn{t}{t} statistics obtained from each of the two
  variables with a given eigenfunction. The \eqn{\tau}{tau} statistic
  can take both positive and negative values, thereby allowing one to
  perform one-directional or two-directional testing.
}
\details{
  \tabular{ll}{
    Package: \tab codep\cr
    Type: \tab Package\cr
    Version: \tab 0.2-1\cr
    Date: \tab 2013-08-06\cr
    License: \tab Copyleft\cr
    LazyLoad: \tab yes\cr
  }
  Function \code{\link{MCA}} performs Multiscale Codependence Analysis
  (MCA).

  Functions \code{\link{test.mca}} and \code{\link{permute.mca}} handle
  parametric permutational testing of the codependence coefficients,
  respectively.

  Methods are provided to print and plot \code{\link{mca-class}} objects
  (\code{\link{print.mca}} and \code{\link{plot.mca}}, respectively) as
  well as summary (\code{\link{summary.mca}}), fitted values
  (\code{\link{fitted.mca}}), residuals (\code{\link{residuals.mca}}),
  and to make predictions (\code{\link{predict.mca}}).

  Function \code{\link{eigenmap}} calculates spatial eigenvector maps
  following the approach outlined in Dray et al. (2006), and which are
  necessary to calculate \code{\link{MCA}}. It returns a
  \code{\link{eigenmap-class}} object. The package also features methods
  to print (\code{\link{print.eigenmap}}) and plot
  (\code{\link{plot.eigenmap}}) these objects.
  
  The package also features an examplary dataset \code{\link{Salmon}}
  containing  76 sampling site positions along a 1520 m river segment as
  well as functions \code{\link{cthreshold}} and
  \code{\link{minpermute}}, which calculates the testwise type I error rate threshold
  corresponding to a give familywise threshold and the number of
  permutations suitable for testing Multiscale Codependence Analysis,
  respectively.
}
\references{
  Dray, S.; Legendre, P. and Peres-Neto, P. 2006. Spatial modelling: a
  comprehensive framework for principal coordinate analysis of neighbor
  matrices (PCNM). Ecol. Modelling 196: 483-493
  
  Guénard, G., Legendre, P., Boisclair, D., and Bilodeau, M. 2010.
  Multiscale codependence analysis: an integrated approach to analyse
  relationships across scales. Ecology 91: 2952-2964
}
\seealso{
  Legendre, P. and Legendre, L. 2012. Numerical Ecology, 3rd English
  edition. Elsevier Science B.V., Amsterdam, The Neatherlands.
}
\author{ Guillaume Guénard, Département des sciences biologiques,
  Université de Montréal, Montréal, Québec, Canada. }
