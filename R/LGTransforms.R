## **************************************************************************
##
##    (c) 2018-2023 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    **Legendre and Gallagher transformation function**
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
#' Transformation for Species Abundance Data
#' 
#' Calculates the transformed species abundances following Legendre and
#' Gallagher.
#' 
#' @param x A species abundance matrix (rows: sites, columns: species).
#' @param method The transformation method, one of "chord" (the default),
#' "chisq", "profile", or "Hellinger" (see details).
#' @param offset Offset value applied to all the columns of \code{x} prior to
#' the other transformations (default: \code{0}, see Details).
#' @param power Exponent for the power transformation (Box-Cox) applied
#' to all columns of \code{x} after the offset and before the transformation
#' specified by argument \code{method} (default: \code{1}, see Details).
#' 
#' @return A matrix of the transformed species abundances.
#' 
#' @details These transformations of species abundances values are useful for
#' multivariate least squares methods for ordination methods, such as the
#' principal component analysis, or modelling methods, such as the multiscale
#' codependence analysis (\code{\link{MCA}}), the canonical redundancy analysis
#' (RDA). They allow one to use least squares methods, which operate on the
#' basis of the Euclidean metric, on species abundance data, for which the
#' Euclidean metric have generally inadequate properties (see Legendre &
#' Gallagher 2001 and Legendre & Borcard 2018, in references below, for a
#' thorough discussion on the topic).
#' 
#' The power (Box Cox) transformation involves the following equation:
#' 
#' y' = (y + offset)^power    if power != 0
#' 
#' y' = log(y + offset)       if power == 0
#' 
#' The default values for \code{offset} (\code{0}) and \code{power} (\code{1})
#' correspond to applying no transformation besides that specified by argument
#' \code{methods}.
#' 
#' @author \packageAuthor{codep}
#' Maintainer: \packageMaintainer{codep}
#' 
#' @references
#' Legendre P. & Gallagher E. D. 2001. Ecologically meaningful transformations
#' for ordination of species data. Oecologia 129: 271-280
#' doi: 10.1007/s004420100716
#' 
#' Box G. E. P. & Cox D. R. 1964. An analysis of transformations. Journal of
#' the Royal Statistical Society Series B 26: 211-243
#' 
#' Legendre P. & Borcard D. 2018. Box-Cox-chord transformations for community
#' composition data prior to beta diversity analysis. Ecography 41: 1820-1824.
#' doi: 0.1111/ecog.03498
#' 
#' Legendre, P. & Legendre, L. 2012. Numerical Ecology, Third English Edition.
#' Elsevier B. V. Amsterdam, The Netherlands.
#' 
#' @importFrom stats dist
#' @importFrom graphics mtext
#' 
#' @examples
#' 
#' data(Doubs)
#' 
#' ## Removing any species that have not been not observed:
#' Doubs.fish -> x
#' x[rowSums(x)!=0,] -> x
#' 
#' ## Transforming the abundances
#' LGTransforms(x,"chord") -> chord
#' LGTransforms(x,"chord",offset=1,power=0) -> log.chord
#' LGTransforms(x,"chord",power=0.25) -> pow.chord
#' LGTransforms(x,"chisq") -> chisq
#' LGTransforms(x,"profile") -> sp_pr
#' LGTransforms(x,"Hellinger") -> Helli
#' 
#' dist(chord)
#' dist(log.chord)
#' dist(pow.chord)
#' dist(chisq)
#' dist(sp_pr)
#' dist(Helli)
#' 
#' ## Legendre & Gallagher synthetic examples:
#' 
#' data(LGDat)
#' 
#' ## Diastemograms:
#' 
#' as.matrix(dist(LGDat[,1L])) -> geo
#' geo[upper.tri(geo)] -> geo
#' 
#' ## Raw Euclidean distances
#' par(mfrow=c(1,1), mar=c(5,5,4,2))
#' 
#' as.matrix(dist(LGDat[,-1L])) -> eco
#' eco[upper.tri(eco)] -> eco
#' 
#' plot(eco~geo, data=data.frame(geo=geo, eco=eco),
#'      xaxp=c(1,18,17), las=1, ylim=c(0,max(eco)),
#'      xlab="True geographic distance",
#'      ylab="Euclidean distance")
#' 
#' ## Euclidean distances on the transformed abundances:
#' par(mfrow=c(3,2), mar=c(3,5,4,2))
#' 
#' LGTransforms(LGDat[,-1L],"chord") -> chord
#' as.matrix(dist(chord)) -> eco
#' eco[upper.tri(eco)] -> eco
#' plot(eco~geo,data=data.frame(geo=geo,eco=eco),
#'      xaxp=c(1,18,17),las=1,xlab="",ylab="",
#'      main="Chord distance",ylim=c(0,max(eco)))
#' 
#' LGTransforms(LGDat[,-1L],"chord",offset=1,power=0) -> log.chord
#' as.matrix(dist(log.chord)) -> eco
#' eco[upper.tri(eco)] -> eco
#' plot(eco~geo,data=data.frame(geo=geo,eco=eco),
#'      xaxp=c(1,18,17),las=1,xlab="",ylab="",
#'      main="Chord distance (log(x+1))",ylim=c(0,max(eco)))
#' 
#' par(mar=c(4,5,3,2))
#' 
#' LGTransforms(LGDat[,-1L],"chord",power=0.25) -> pow.chord
#' as.matrix(dist(pow.chord)) -> eco
#' eco[upper.tri(eco)] -> eco
#' plot(eco~geo,data=data.frame(geo=geo,eco=eco),
#'      xaxp=c(1,18,17),las=1,xlab="",ylab="",
#'      main="Chord distance (power=0.25)",ylim=c(0,max(eco)))
#' 
#' LGTransforms(LGDat[,-1L],"chisq") -> chisq
#' as.matrix(dist(chisq)) -> eco
#' eco[upper.tri(eco)] -> eco
#' plot(eco~geo,data=data.frame(geo=geo,eco=eco),
#'      xaxp=c(1,18,17),las=1,xlab="",ylab="",
#'      main="Chi-square distance",ylim=c(0,max(eco)))
#' 
#' par(mar=c(5,5,2,2))
#' 
#' LGTransforms(LGDat[,-1L],"profile") -> sp_pr
#' as.matrix(dist(sp_pr)) -> eco
#' eco[upper.tri(eco)] -> eco
#' plot(eco~geo,data=data.frame(geo=geo,eco=eco),
#'      xaxp=c(1,18,17),las=1,xlab="",ylab="",
#'      main="Dist. between profiles",ylim=c(0,max(eco)))
#' 
#' LGTransforms(LGDat[,-1L],"Hellinger") -> Helli
#' as.matrix(dist(Helli)) -> eco
#' eco[upper.tri(eco)] -> eco
#' plot(eco~geo,data=data.frame(geo=geo,eco=eco),
#'      xaxp=c(1,18,17),las=1,xlab="",ylab="",
#'      main="Hellinger distance",ylim=c(0,max(eco)))
#' 
#' mtext(text="True geographic distance", side=1, line=-1.5, outer=TRUE)
#' mtext(text="Ecological distance", side=2, line=-1.5, outer=TRUE)
#' 
#' ## Examples from Legendre & Legendre 2012, page 329 (Figure 7.8):
#' 
#' matrix(c(0,0,1,4,1,0,8,1,0),3L,3L) -> LL329
#' 
#' ## D1:  Euclidean distance
#' dist(LL329)
#' 
#' ## Chord transformation (D3: chord distance)
#' LGTransforms(LL329,"chord") -> tr       
#' tr
#' dist(tr)
#' 
#' ## "Species profile" transformation (D18)
#' LGTransforms(LL329,"profile") -> tr       
#' tr
#' dist(tr)
#' 
#' ## Hellinger transformation (D17: Hellinger distance)
#' LGTransforms(LL329,"Hellinger") -> tr       
#' tr
#' dist(tr)
#' 
#' ## Chi-square transformation (D16: Chi-square distance)
#' LGTransforms(LL329,"chisq") -> tr       
#' tr
#' dist(tr)
#' 
#' @useDynLib codep, .registration = TRUE
#' 
#' @export
LGTransforms <- function(x, method = c("chord","chisq","profile","Hellinger"),
                         offset = 0, power = 1) {
  method <- match.arg(method)
  method <- match(method, c("chord","chisq","profile","Hellinger"))
  if(!is.matrix(x))
    x <- as.matrix(x)
  if(offset != 0)
    x <- x + offset
  if(power != 1) {
    if(power != 0) {
      x <- x^power
    } else {
      x <- log(x)
    }
  }
  storage.mode(x) <- "double"
  return(
    .C(
      "LGTr_C",
      x,
      nrow(x),
      ncol(x),
      method
    )[[1L]]
  )
}
#' 
