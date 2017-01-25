\encoding{utf8}
\name{Mite}
\docType{data}
\alias{Mite}
\alias{mite.species}
\alias{mite.env}
\alias{mite.geo}
\title{Lac Geai oribatid mites community data}
\description{
  A data set containing information on the microfauna of 70 sites
  located within the the peat on shore of Lac Geai, a bog lake located
  at the Station de Biologie de l'Université de Montréal, St-Hippolyte,
  QC, Canada (+45.9954; -73.9936).
}
\usage{Mite}
\format{ Contains three matrices:
  \code{mite.species} the abundance of 35 morpho-species of oribatid
  mites (Acari).

  \code{mite.env} 14 environmental variables (quantitative and binary).

  \code{mite.geo} the relative coordinates of the samples.
}
\details{
  Values in \code{mite.species} are counts of individuals of each of the
  morpho-species obtained from 5 cm diameter cores going from the
  surface of the peat down to a depth of 7 cm. See Bordard & Legendre
  (1994) and reference therein for details about sample treatment and
  species identification.

  \code{mite.env} contains two quantitative variables, namely the
  substratum density (g/L) and water content (percent wet mass over dry
  mass), in addition to 12 dummy variables. The first seven represent
  the composition of the substratum: \emph{Sphagnum magellacinum} (with
  a majority of \emph{S. rubellum}), \emph{S. rubellum},
  \emph{S. nemorum}, (with a majority of \emph{S. augustifollium}),
  \emph{S. rubellum} + \emph{S. magellicum} (in equal proportions),
  lignous litter, bare peat, and interface between \emph{Sphagnum}
  species. The next three dummy variables represent the presence and
  abundance of shrubs (\emph{Kalmia polifolia}, \emph{K. angustifolia},
  and \emph{Rhododentron groenlandicum}): none, few, and many. The last
  two dummy variables represent the microtopography of the peat: blanket
  (flat) or hummock (raised).

  \code{mite.geo} contains the location of the samples, in meters, with
  respect to the sampling grid. Point (0,0) is the lower left end of the
  plot for an observer looking from the shore towards the water. The
  \code{x} coordinate is the offset along the shore (from left to right)
  while the \code{y} coordinate is the offset from the shore while
  moving towards the water (See Borcard & Legendre, 1994, Fig. 1 for
  details on the sampling area).
}
\source{ Daniel Borcard, Département de sciences biologiques, Université
  de Montréal, Montréal, Québec, Canada. }
\references{
  Borcard, D. & Legendre, P. 1994. Environmental control and spatial
  structure in ecological communities: an example using Oribatid mites
  (Acari, Oribatei). Environ. Ecol. Stat. 1: 37-61
}
\seealso{
  Borcard, D.; P. Legendre & P. Drapeau. 1992. Partialling out the
  spatial component of ecological variation. Ecology 73: 1045-1055

  Legendre, P. 2005. Species associations: the Kendall coefficient of
  concordance revisited. Journal of Agricultural, Biological and
  Environmental Statistics 10: 226-245

  Borcard, D.; Gillet, F. & Legendre, P. 2011. Numerical Ecology with
  R. Springer, New-York, NY, USA.
}
\examples{
data(Mite)
summary(mite.species)
summary(mite.env)
summary(mite.geo)
}
\keyword{ Mite }
