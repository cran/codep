## **************************************************************************
##
##    (c) 2018-2023 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    **Legendre and Gallagher synthetic example data set**
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
##    Data set documentation generation file
##
## **************************************************************************
##
#' Legendre and Gallagher Synthetic Example
#' 
#' A data set used as a synthetic example in paper Legendre and Gallagher
#' (2001).
#' 
#' @docType data
#' 
#' @name LGDat
#' 
#' @usage data(LGDat)
#' 
#' @format A 19 rows by 10 columns \code{\link{data.frame}}.
#' 
#' @details This synthetic data set is described by Legendre and Gallagher
#' (2001) and was used to test species abundance transformations. Its first
#' column contains geographic locations from 1 to 19 (no particular units are
#' specified). The five columns that follow contain abundances of five species
#' with abundances peaking at 7-8 at different locations (site 1, 5, 10, 15, and
#' 19). The latter are considered "abundant species". For next four columns
#' contains abundances of "rare species" occurring in between the abundance
#' species (abundances from 1 to 4).
#' 
#' @source Legendre, P. & Gallagher E. D. 2001. Ecologically meaningful
#' transformations for ordination of species data. Oecologia 129: 271-280
#' doi: 10.1007/s004420100716
#' 
#' @examples
#' data(LGDat)
#' summary(LGDat)
#' 
NULL
#' 
