## **************************************************************************
##
##    (c) 2018-2023 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    **Oribatid mite data set**
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
#' The Oribatid Mite Data Set
#' 
#' Borcard et al's oribatid mite community composition from Lac Geai, Canada.
#' 
#' @docType data
#' 
#' @name mite
#' 
#' @aliases mite.species mite.env mite.geo
#' 
#' @usage data(mite)
#' 
#' @format Contains three matrices:
#' \describe{
#' \item{mite.species}{ The abundance of 35 morpho-species of oribatid mites
#' (Acari). }
#' \item{mite.env}{ 14 environmental variables (quantitative and binary). }
#' \item{mite.geo}{ The relative coordinates of the samples. }
#' }
#' 
#' @details Values in \code{mite.species} are counts of individuals of each of
#' the morpho-species obtained from 5 cm diameter cores going from the surface
#' of the peat down to a depth of 7 cm. See Bordard & Legendre (1994) and
#' reference therein for details about sample treatment and species
#' identification.
#' 
#' `mite.env` contains two quantitative variables, namely the substratum density
#' (g/L) and water content (percent wet mass over dry mass), in addition to 12
#' dummy variables. The first seven represent the composition of the substratum:
#' \emph{Sphagnum magellacinum} (with a majority of \emph{S. rubellum}),
#' \emph{S. rubellum}, \emph{S. nemorum}, (with a majority of
#' \emph{S. augustifollium}), \emph{S. rubellum} + \emph{S. magellicum} (in
#' equal proportions), lignous litter, bare peat, and interface between
#' \emph{Sphagnum} species. The next three dummy variables represent the
#' presence and abundance of shrubs (\emph{Kalmia polifolia},
#' \emph{K. angustifolia}, and \emph{Rhododentron groenlandicum}): none, few,
#' and many. The last two dummy variables represent the microtopography of the
#' peat: blanket (flat) or hummock (raised).
#' 
#' `mite.geo` contains the location of the samples, in meters, with respect to
#' the sampling grid. Point (0,0) is the lower left end of the plot for an
#' observer looking from the shore towards the water. The `x` coordinate is the
#' offset along the shore (from left to right) while the `y` coordinate is the
#' offset from the shore while moving towards the water (See Borcard &
#' Legendre, 1994, Fig. 1 for details on the sampling area).
#' 
#' @source Daniel Borcard, Département de sciences biologiques, Université de
#' Montréal, Montréal, Québec, Canada.
#' 
#' @references
#' Borcard, D. & Legendre, P. 1994. Environmental control and spatial structure
#' in ecological communities: an example using Oribatid mites (Acari, Oribatei).
#' Environ. Ecol. Stat. 1: 37-61
#' 
#' @seealso
#' Borcard, D.; P. Legendre & P. Drapeau. 1992. Partialling out the spatial
#' component of ecological variation. Ecology 73: 1045-1055
#' 
#' Legendre, P. 2005. Species associations: the Kendall coefficient of
#' concordance revisited. Journal of Agricultural, Biological and Environmental
#' Statistics 10: 226-245
#' 
#' Borcard, D.; Gillet, F. & Legendre, P. 2011. Numerical Ecology with R.
#' Springer, New-York, NY, USA.
#' 
#' @examples
#' data(mite)
#' summary(mite.species)
#' summary(mite.env)
#' summary(mite.geo)
#' 
#' @keywords mite
NULL
#' 
