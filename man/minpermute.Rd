\encoding{utf8}
\name{minpermute}
\alias{minpermute}
\title{ Number of permutations for MCA }
\description{
  Calculate the number of permutations suitable for testing Multiscale
  Codependence Analysis.
}
\usage{
minpermute(alpha,nbtest,margin=1,ru=3)
}
\arguments{
  \item{alpha}{ The familywise type I error threshold allowable for the
    complete analysis. }
  \item{nbtest}{ The number of test performed (the number of
    eigenvectors in the \sQuote{mem} object in the case of
    \code{\link{MCA}}). }
  \item{margin}{ A margin allowed for the number of permutation. Default
    value: 1. }
  \item{ru}{ The magnitude of the round-up to apply to the number of
    permutations. }
}
\value{
  The minimum number of permutation to be used for
  \code{\link{permute.cdp}}.
}
\details{
  This function calculate the number of permutations for use with
  \code{\link{permute.cdp}}. Parameter \code{margin} allows to apply a
  safe margin to the number of permutations. The minimal suitable value
  for this parameter is 1. Parameter \code{ru} allows one to round-up
  the number of permutations. A value of 0 implies no round-up, a value
  of 1 a round-up to the next ten, 2 a round-up to the next hundred, and
  so on. Function \code{minpermute} is called internally by
  \code{\link{permute.cdp}} in case \code{permute = NA}. In that case,
  the margin is set to 10 (\code{margin = 10}) and the outcome is
  rounded-up to the next thousand (\code{ru = 3}). This function is
  meant for users that wish to apply their own margins and round-up
  factors to calculate the number of permutations for use with
  \code{permute.cdp}.
}
\references{
  Guénard, G., Legendre, P., Boisclair, D., and Bilodeau, M. 2010.
  Multiscale codependence analysis: an integrated approach to analyse
  relationships across scales. Ecology 91: 2952-2964
}
\seealso{
  \link{permute.cdp}
}
\examples{
# For a 5\% threshold under 50 tests.
minpermute(alpha = 0.05, nbtest=50)
# Allowing more margin (implies more computation time).
minpermute(alpha = 0.05, nbtest=50, margin=10, ru=3)
}
\author{ Guillaume Guénard, Département des sciences biologiques,
  Université de Montréal, Montréal, Québec, Canada. }
