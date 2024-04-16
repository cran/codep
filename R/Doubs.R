## **************************************************************************
##
##    (c) 2018-2023 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    **Doubs data set**
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
#' The Doubs Fish Data
#' 
#' Fish community composition of the Doubs River, France.
#' 
#' @docType data
#' 
#' @name Doubs
#' 
#' @aliases Doubs.fish Doubs.env Doubs.geo
#' 
#' @usage data(Doubs)
#' 
#' @format Contains three matrices:
#' \describe{
#' \item{Doubs.fish}{ The abundance of 27 fish species. }
#' \item{Doubs.env}{ Nine environmental variables (all quantitative). }
#' \item{Doubs.geo}{ Geographic information of the samples. }
#' }
#' 
#' @details Values in `Doubs.fish` are counts of individuals of each of 27
#' species observed in a set of 30 sites located along the 453 km long Doubs
#' River, France (see Verneaux 1973 for further details about fishing methods
#' and effort).
#' \describe{
#'   \item{Doubs.env}{ contains 11 quantitative variables, namely the slope
#'   (`slo`; 1/1000) and mean minimum discharge (`flo` m³/s) of the river, the
#'   pH of the water, its harness (Calcium concentration; `har`; mg/L),
#'   phosphate (`pho`; mg/L), nitrate (`nit`; mg/L), and ammonium (`amm`; mg/L),
#'   concentration as well as its dissolved oxygen (`oxy`; mg/L) and biological
#'   oxygen demand (`bdo`; mg/L).}
#' \item{Doubs.geo}{ contains geographical information. `Lon`, the longitude and
#'   `Lat`, the latitude of the sample (degree) as well as `DFS`, its distance
#'   from the source of the river (km) and `Alt`, altitude (m above see
#'   level). }
#' }
#' 
#' @source Verneaux, 1973
#' 
#' @references
#' Verneaux J. 1973. - Cours d'eau de Franche-Comté (Massif du Jura). Recherches
#' écologiques sur le réseau hydrographique du Doubs. Essai de biotypologie.
#' Thèse d'état, Besançon. 257 p.)
#' 
#' Verneaux, J.; Schmitt, V.; Verneaux, V. & Prouteau, C. 2003. Benthic insects
#' and fish of the Doubs River system: typological traits and the development of
#' a species continuum in a theoretically extrapolated watercourse.
#' Hydrobiologia 490: 60-74
#' 
#' @seealso
#' Borcard, D.; Gillet, F. & Legendre, P. 2011. Numerical Ecology with R.
#' Springer, New-York, NY, USA.
#' 
#' @examples
#' data(Doubs)
#' summary(Doubs.fish)
#' summary(Doubs.env)
#' summary(Doubs.geo)
#' 
#' @keywords Doubs
NULL
#' 
