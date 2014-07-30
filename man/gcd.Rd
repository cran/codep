\encoding{utf8}
\name{gcd}
\alias{gcd.slc}
\alias{gcd.hf}
\alias{gcd.vife}
\title{Great circle distances}
\description{
  Functions to calculate great circle distances among sets of points
  with known longitude and lattitudes.
}
\usage{
gcd.slc(ss, radius = 6371)
gcd.hf(ss, radius = 6371)
gcd.vife(ss, a = 6378137, b = 6356752.314245, f = 1/298.257223563)
}
\arguments{
  \item{ss}{ A two-columns data.frame or matrix of the geographic
    coordinates (in degree.decimals; lattitude and then longitude) of
    the locations between which the geodesic distances are being
    calculated. }
  \item{radius}{ Mean earth radius (mean length of parallels) in
    km. Default: 6371 km }
  \item{a}{ Length (in m) of major axis of the ellipsoid (radius at
    equator). Default: 6378137 m (i.e. that for WGS-84). }
  \item{b}{ Length (in m) of minor axis of the ellipsoid (radius at the
    poles). Default: 6356752.314245 m (i.e. that for WGS-84). }
  \item{f}{ Flattening of the ellipsoid. Default: 1/298.257223563
    (i.e. that for WGS-84). }
}
\details{
 The calculation of spatial eigenvector maps requires a distance
  matrix. The euclidean distance is appropriate when a cartesian plan
  can be reasonably assumed. The latter can be calculated with function
  \code{\link{dist}}. The great circle distance is appropriate when the
  sampling points are located on a spheroid (e.g. planet
  Earth). Function \code{gcd.slc} uses the spherical law of cosines,
  which is the fastest of the three approaches and performs relatively
  well for distances above 1 m. Function \code{gcd.hf} uses the the
  Haversine formula, which is more accurate for smaller distances
  (bellow 1 m). Finally, functions \code{gcd.vife} uses the Vincenty
  inverse formula for ellipsoids, which is an iterative approach that
  take substantially more computation time than the latter two methods
  has precision up to 0.5 mm with exact longitudes and lattitudes. }
\value{
  An object of class ‘"dist"’ (see \code{\link{dist}} for details) that
  contains the distances (in km) between the locations (rows of
  \code{ss}).
}
\references{
  Mario Pineda-Krch, URL:
  http://www.r-bloggers.com/great-circle-distance-calculations-in-r/

  Vincenty, T. (1975) Closed formulas for the direct and reverse
    geodetic problems. J. Geodesy 51(3): 241-342
}
\author{ Guillaume Guénard, Departement des sciences biologiques,
  Universite de Montréal, Montréal, Quebec, Canada. }
\seealso{
  \code{\link{eigenmap-class}}
}
\examples{
#
# Calculating the distances between Canada's capital cities:
CapitalCitiesOfCanada <-
    matrix(c(45.417,-75.7,53.533333,-113.5,48.422151,-123.3657,
             49.899444,-97.139167,45.95,-66.666667,47.5675,-52.707222,
             44.647778,-63.571389,43.7,-79.4,46.24,-63.1399,
             46.816667,-71.216667,50.454722,-104.606667,62.442222,-114.3975,
             63.748611,-68.519722,60.716667,-135.05),14L,2L,byrow=TRUE,
             dimnames=list(c("Ottawa","Edmonton","Victoria","Winnipeg",
                             "Fredericton","St-John's","Halifax","Toronto",
                             "Charlottetown","Quebec City","Regina",
                             "Yellowknife","Iqaluit","Whitehorse"),c("Lon","Lat")))
#
sphericalcosdists <- gcd.slc(CapitalCitiesOfCanada)
vincentydists <- gcd.vife(CapitalCitiesOfCanada)
#
cor(as.numeric(sphericalcosdists),as.numeric(vincentydists))
percentdev <- 100*(vincentydists-sphericalcosdists)/vincentydists
mean(percentdev)
# Spherical Law of Cosines underestimated these distances by ~0.26
# percent.
#
}
