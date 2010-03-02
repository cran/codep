% --- Source file: ./man/cthreshold.Rd ---
% File src/library/codep/man/cthreshold.Rd
% Part of the R package, http://www.R-project.org
% Copyright 2009 Guenard, Guillaume and Legendre, Pierre
% Distributed under GPL 2 or later
\name{cthreshold}
\alias{cthreshold}
\title{ Familywise type I error rate }
\description{
  Function to calculate the testwise type I error rate threshold that
  corresponds to a give familywise threshold.
}
\usage{
cthreshold(alpha, nbtest)
}
\arguments{
  \item{alpha}{ The familywise type I error threshold. }
  \item{nbtest}{ The number of tests performed. }
}
\value{
  The threshold that have to be used for individual tests.
}
\details{
  Type I error rate inflation occurs when a single hypothesis is tested
  indirectly using inferences about two or more (\emph{i.e.} a family
  of) sub-hypotheses. In such situation, the probability of type I error
  (\emph{i.e.} the probability of incorrectly rejecting the null
  hypothesis) of the single, familywise, hypothesis is higher than the
  lowest, testwise, probabilities. As a consequence, the rejection of
  null hypothesis for one or more individual tests does not warrant that
  the correct decision (whether to reject the the null hypothesis on a
  familywise basis) was taken properly. This function allows to obtain
  correct, familywise, alpha thresholds in the context of multiple
  testing. It is base on the Sidak inegality.
}
\references{
  Sidak, Z. 1967. Rectangular Confidence Regions for Means of
  Multivariate Normal Distributions J. Am. Stat. Assoc. 62: 626-633
  
  Wright, P. S. 1992. Adjusted p-values for simultaneous inference.
  Biometrics 48: 1005-1013  
}
\seealso{
  Legendre, P. and Legendre, L. 1998. Numerical Ecology. Elsevier
  Science B.V., Amsterdam, The Neatherlands. p. 18
}
\examples{
# For a familywise threshold of 5\% with 5 tests:
cthreshold(c(0.05),5)   # Threshold of each test is 0.01020622
}
\author{ Guillaume Guenard, Laboratoire evolution et diversite
  biologique, CNRS / Universite Paul-Sabatier, Toulouse, France. }