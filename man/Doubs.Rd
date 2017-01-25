\encoding{utf8}
\name{Doubs}
\docType{data}
\alias{Doubs}
\alias{Doubs.fish}
\alias{Doubs.env}
\alias{Doubs.geo}
\title{The Doubs fish data}
\description{
  Fish community composition of the Doubs River, France.
}
\usage{Doubs}
\format{ Contains three matrices:
  \code{Doubs.fish} the abundance of 27 fish species.

  \code{Doubs.env} 9 environmental variables (all quantitative).

  \code{Doubs.geo} geographic information of the samples.
}
\details{
  Values in \code{Doubs.fish} are counts of individuals of each of 27
  species observed in a set of 30 sites located along the 453 km long
  Doubs River, France (see Verneaux 1973 for further details about
  fishing methods and effort).
  
  \code{Doubs.env} contains 11 quantitative variables, namely the slope
  (\code{slo} 1/1000) and mean minimum discharge (\code{flo} m³/s) of
  the river, the pH of the water, its harness (Calcium concentration;
  \code{har}; mg/L), phosphate (\code{pho}; mg/L), nitrate (\code{nit};
  mg/L), and ammonium (\code{amm}; mg/L), concentration as well as its
  dissolved oxygen (\code{oxy}; mg/L) and biological oxygen demand
  (\code{bdo}; mg/L).
  
  \code{Doubs.geo} contains geographical information. \code{Lon}, the
  longitude and \code{Lat}, the latitude of the sample (degree) as well
  as \code{DFS}, its distance from the source of the river (km) and
  \code{Alt}, altitude (m above see level).
}
\source{ Verneaux, 1973 }
\references{
  Verneaux J. 1973. - Cours d'eau de Franche-Comté (Massif du
  Jura). Recherches écologiques sur le réseau hydrographique du
  Doubs. Essai de biotypologie. Thèse d'état, Besançon. 257 p.)

  Verneaux, J.; Schmitt, V.; Verneaux, V. & Prouteau, C. 2003. Benthic
  insects and fish of the Doubs River system: typological traits and the
  development of a species continuum in a theoretically extrapolated
  watercourse. Hydrobiologia 490: 60-74
}
\seealso{
  Borcard, D.; Gillet, F. & Legendre, P. 2011. Numerical Ecology with
  R. Springer, New-York, NY, USA.
}
\examples{
data(Doubs)
summary(Doubs.fish)
summary(Doubs.env)
summary(Doubs.geo)
}
\keyword{ Doubs }
