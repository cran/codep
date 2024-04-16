## **************************************************************************
##
##    (c) 2018-2023 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    **cthreshold function**
##
##    This file is part of codep
##
##    codep is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    codep is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with codep. If not, see <https://www.gnu.org/licenses/>.
##
##    R source code file
##
## **************************************************************************
##
#' Familywise Type I Error Rate
#' 
#' Function to calculate the testwise type I error rate threshold corresponding
#' to a give familywise threshold.
#' 
#' @param alpha The familywise type I error threshold.
#' @param nbtest The number of tests performed.
#' 
#' @return The threshold that have to be used for individual tests.
#' 
#' @details Type I error rate inflation occurs when a single hypothesis is
#' tested indirectly using inferences about two or more (i.e., a family of)
#' sub-hypotheses. In such situation, the probability of type I error (i.e., the
#' probability of incorrectly rejecting the null hypothesis) of the single,
#' familywise, hypothesis is higher than the lowest, testwise, probabilities. As
#' a consequence, the rejection of null hypothesis for one or more individual
#' tests does not warrant that the correct decision (whether to reject the the
#' null hypothesis on a familywise basis) was taken properly. This function
#' allows to obtain correct, familywise, alpha thresholds in the context of
#' multiple testing. It is base on the Sidak inegality.
#' 
#' @author \packageAuthor{codep}
#' Maintainer: \packageMaintainer{codep}
#' 
#' @seealso Legendre, P. and Legendre, L. 1998. Numerical Ecology. Elsevier
#' Science B.V., Amsterdam, The Neatherlands. p. 18
#' 
#' @references
#' Sidak, Z. 1967. Rectangular Confidence Regions for Means of Multivariate
#' Normal Distributions J. Am. Stat. Assoc. 62: 626-633
#' 
#' Wright, P. S. 1992. Adjusted p-values for simultaneous inference. Biometrics
#' 48: 1005-1013
#' 
#' @examples
#' ## For a familywise threshold of 5% with 5 tests:
#' cthreshold(0.05, 5)   ## The corrected threshold for each test is 0.01020622
#' 
#' @export
cthreshold <- function(alpha, nbtest) return(1 - (1 - alpha)^(nbtest^-1))
#' 
