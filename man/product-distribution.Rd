\encoding{utf8}
\name{Product-distribution}
\alias{dtau}
\alias{ptau}
\title{Frequency distributions for MCA parametric testing}
\description{
  Density and distribution functions of the tau statistic, which is the
  product of two Student's t distributions with the same number of
  degrees of freedom.
}
\usage{
dtau(x, nu, tol = .Machine$double.eps ^ 0.5)
ptau(q, nu, lower.tail = TRUE, tol = .Machine$double.eps ^ 0.5)
}
\arguments{
  \item{x, q}{ vector of quantile. }
  \item{nu}{ degrees of freedom (>0, maybe non-integer). \code{nu = Inf}
    is allowed. }
  \item{lower.tail}{ logical; if TRUE (default), probabilities are P[X
    <= x], otherwise, P[X > x]. }
  \item{tol}{ tolerance used for numerical estimation. }
}
\value{
  \code{dtau} gives the density and \code{ptau} gives the distribution
  function.
}
\details{
  The density distribution of a variable \code{z} that is the product of
  two random variables \code{x} and \code{y} with density distributions
  f(x) and g(y) is the integral over the intersection of the domains of
  \code{x} and \code{y} of f(x) * g(z/x) / abs(x) dx.

  \code{dtau} estimate density values using numerical integration
  (\code{\link{integrate}}) and Student's t density distribution
  function \code{\link{dt}} with 'df = nu'. \code{ptau} integrates
  \code{dtau} on the interval [\code{-Inf},\code{q}] when 'lower.tail =
  TRUE' (the default) and on the interval [\code{q},\code{Inf}] when
  'lower.tail = FALSE'.
}
\references{
  Springer, M. D. 1979. The algebra of random variables. John Wiley and
  Sons Inc., Hoboken, NJ, USA.

  Guénard, G., Legendre, P., Boisclair, D., and Bilodeau, M. 2010.
  Multiscale codependence analysis: an integrated approach to analyse
  relationships across scales. Ecology 91: 2952-2964
}
\seealso{
  \link{test.mca}
}
\examples{
#
### Displays the probability distribution function for five different
### numbers of degrees of freedom:
#
x <- c(-(10^seq(0.5,-4,-0.05)),10^seq(-4,0.5,0.05))
plot(y=dtau(x,1),x=x,type="l",col="black",las=1,ylab="Density",ylim=c(0,0.5))
lines(y=dtau(x,2),x=x,col="purple")
lines(y=dtau(x,5),x=x,col="blue")
lines(y=dtau(x,10),x=x,col="green")
lines(y=dtau(x,100),x=x,col="red")
#
### Displays the density distribution function for 10 degrees of freedom
### and the cumulative probability below x = -1 and above x = 1.
#
y <- dtau(x,10)
plot(y=y,x=x,type="l",col="black",las=1,ylab="Density",ylim=c(0,0.5))
polygon(x=c(x[which(x==1):length(x)],x[length(x)],1),
        y=c(y[which(x==1):length(x)],0,0),col="grey")
text(round(ptau(1,10,lower.tail=FALSE),3),x=1.5,y=0.03)
polygon(x=c(-1,x[1L],x[1L:which(x==-1)]),y=c(0,0,y[1L:which(x==-1)]),col="grey")
text(round(ptau(-1,10),3),x=-1.5,y=0.03)
#
}
\author{ Guillaume Guénard, Département des sciences biologiques,
  Université de Montréal, Montréal, Québec, Canada. }
