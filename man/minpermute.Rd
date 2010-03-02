% --- Source file: ./man/minpermute.Rd ---
% File src/library/codep/man/minpermute.Rd
% Part of the R package, http://www.R-project.org
% Copyright 2009 Guenard, Guillaume and Legendre, Pierre
% Distributed under GPL 2 or later
\name{minpermute}
\alias{minpermute}
\title{ Number of permutations for MCA }
\description{
  Calculate the number of permutations suitable for Multiscale
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
    \code{\link{mca}}). }
  \item{margin}{ A margin allowed for the number of permutation. }
  \item{ru}{ The magnitude of the round-up to apply to the number of
    permutations. }
}
\value{
  The minimum number of permutation to be used for
  \code{\link{permute.mca}}.
}
\details{
  This function calculate the number of permutations for use with
  \code{\link{permute.mca}}. Parameter \code{margin} allows to apply a
  safe margin to the number of permutations. The minimal suitable value
  for this parameter is 1. Parameter \code{ru} allows to round-up the
  number of permutations. A value of 0 implies no round-up, a value of 1
  a round-up to the next ten, 2 a round-up to the next hundred, and so
  on. Function \code{minpermute} is called internally by
  \code{\link{permute.mca}} in case \code{permute = NA}. In that case,
  the margin is set to 10 (\code{margin = 10}) and the outcome is
  rounded-up the next thousand (\code{ru = 3}). This function is meant
  for users that wish to apply their own margins and round-up factors to
  calculate the number of permutations for use with \code{permute.mca}.
}
%\references{}
\seealso{
  \link{permute.mca}
}
\examples{
# For a 5\% threshold under 50 tests.
minpermute(alpha = 0.05, nbtest=50)
# Allowing more margin (implies more computation time).
minpermute(alpha = 0.05, nbtest=50, margin=10, ru=3)
}
\author{ Guillaume Guenard, Laboratoire evolution et diversite
  biologique, CNRS / Universite Paul-Sabatier, Toulouse, France. }
