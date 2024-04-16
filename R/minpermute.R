## **************************************************************************
##
##    (c) 2018-2023 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    **minpermute function**
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
#' Number of Permutations for MCA
#' 
#' Calculates the number of permutations suitable for testing Multiscale
#' Codependence Analysis.
#' 
#' @param alpha The familywise type I error threshold allowable for the complete
#' analysis.
#' @param nbtest The number of test performed (the number of eigenvectors in the
#' `mem` object in the case of \code{\link{MCA}}).
#' @param margin A margin allowed for the number of permutation. Default value:
#' 1.
#' @param ru The magnitude of the round-up to apply to the number of
#' permutations.
#' 
#' @return The minimum number of permutation to be used for
#' \code{\link{permute.cdp}}.
#' 
#' @details This function calculate the number of permutations for use with
#' \code{\link{permute.cdp}}. Argument \code{margin} allows to apply a safe
#' margin to the number of permutations. The minimal suitable value for this
#' parameter is 1. Argument \code{ru} allows one to round-up the number of
#' permutations. A value of 0 implies no round-up, a value of 1 a round-up to
#' the next ten, 2 a round-up to the next hundred, and so on. Function
#' \code{\link{minpermute}} is called internally by \code{\link{permute.cdp}} in
#' case \code{permute = NA}. In that case, the margin is set to 10
#' (\code{margin = 10}) and the outcome is rounded-up to the next thousand
#' (\code{ru = 3}). This function is meant for users that wish to apply their
#' own margins and round-up factors to calculate the number of permutations for
#' use with \code{\link{permute.cdp}}.
#' 
#' @author \packageAuthor{codep}
#' Maintainer: \packageMaintainer{codep}
#' 
#' @seealso \link{permute.cdp}
#' 
#' @references
#' Guénard, G., Legendre, P., Boisclair, D., and Bilodeau, M. 2010. Multiscale
#' codependence analysis: an integrated approach to analyse relationships across
#' scales. Ecology 91: 2952-2964
#' 
#' @examples
#' ## For a 5\% threshold under 50 tests.
#' minpermute(alpha = 0.05, nbtest=50)
#' 
#' ## Allowing more margin (implies more computation time).
#' minpermute(alpha = 0.05, nbtest=50, margin=10, ru=3)
#' 
#' @export
minpermute <- function(alpha, nbtest, margin = 1, ru = 3)
  return(
    round(
      floor(
        margin * (1 - (1 - alpha)^(nbtest^-1))^-1
      ),
      -ru
    ) + (10^ru) - 1
  )
#' 
