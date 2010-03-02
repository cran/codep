#
### Computation of Multiscale Codependence Analysis.
### Calculations outlined in Guenard et al. (In prep)
### mca() returns a list object with a 'mca' class attribute.
### Dedicated printing and plotting functions are provided for class 'mca'
### Required package: <none>
###
### Guillaume Guénard, January 2009
#
.First.lib <- function (lib, pkg) {
    library.dynam("codep",pkg, lib)
}

cthreshold <- function(alpha,nbtest) return(1-(1-alpha)^(nbtest^-1))
#
minpermute <- function(alpha,nbtest,margin=1,ru=3) {
  return(round(floor(margin*(1-(1-alpha)^(nbtest^-1))^-1),-ru)+(10^ru)-1)
}
#
mca <- function(y, x, memobj) {
  if(!inherits(memobj,"mem")) stop("Parameter 'memobj' must be a 'mem' object!")
  if(length(y) != length(x)) stop("Number of observations in y and x do not match!")
  if(nrow(memobj$U) != length(y)) stop("Number of observations in y does not match the number of lines in U.")
  #
  mssd <- c(my=NA,mx=NA,ssdy=NA,ssdx=NA)
  mssd[1:2] <- c(mean(y),mean(x))
  yxc <- cbind(yc=y-mssd[1], xc=x-mssd[2])
  mssd[3:4] <- c(t(yxc[,1]) %*% yxc[,1], t(yxc[,2]) %*% yxc[,2])  
  Upyxcb <- matrix(NA, dim(memobj$U)[2], 4)  
  colnames(Upyxcb) <- c("Upy", "Upx", "C", "B")
  rownames(Upyxcb) <- colnames(memobj$U)
  Upyxcb[,1:2] <- t(memobj$U) %*% yxc
  Upyxcb[,3] <- (Upyxcb[,1] * Upyxcb[,2]) / sqrt(mssd[3] * mssd[4])
  Upyxcb[,4] <- (Upyxcb[,1] / Upyxcb[,2])
  #
  ## Output block.
  return(structure(list(
    data = cbind(y=y,x=x),
    memobj = memobj,
    Upyxcb = Upyxcb,
    test = NULL),
    class = "mca"))
}
#
test.mca <- function(mcaobj, alpha = 0.05, max.step = Inf) {
  if(!is.finite(max.step[1])) max.step <- ncol(mcaobj$memobj$U)
  us <- matrix(NA, nrow(mcaobj$memobj$U), 0)
  uspyx <- matrix(NA, 0, 2)
  yxc <- cbind(yc=mcaobj$data[,1]-mean(mcaobj$data[,1]), xc=mcaobj$data[,2]-mean(mcaobj$data[,2]))
  ord <- order(abs(mcaobj$Upyxcb[,3]), decreasing = TRUE)
  ttable <- matrix(NA, 0, 4)
  colnames(ttable) <- c("t2", "ddf", "Testwise p", "Familywise p")
  step <- 1
  while(step != 0) {
    us <- cbind(us, mcaobj$memobj$U[,ord[step]])
    uspyx <- rbind(uspyx, mcaobj$Upyxcb[ord[step],1:2])
    ddfr <- nrow(mcaobj$data) - step - 1
    ryx <- yxc - (us %*% uspyx)
    t2 <- ddfr * (uspyx[step,1] * uspyx[step,2]) / sqrt((t(ryx[,1]) %*% ryx[,1]) * (t(ryx[,2]) %*% ryx[,2]))
    t2s <- ddfr * min(uspyx[step,1]^2 / (t(ryx[,1]) %*% ryx[,1]), uspyx[step,2]^2 / (t(ryx[,2]) %*% ryx[,2]))
    ttable <- rbind(ttable, c(t2, ddfr, NA, NA))
    ttable[step,3] <- pf(t2s,1,ddfr,lower.tail=FALSE)^2
    ttable[step,4] <- 1 - (1 - ttable[step,3])^(ncol(mcaobj$memobj$U) - step + 1)
    if (ttable[step,4] > alpha || step >= max.step) {
      rownames(ttable) <- colnames(mcaobj$memobj$U)[ord[1:step]]
      step <- 0
    }
    else step <- step + 1
  }
  signif <- ord[which(ttable[,4] <= alpha)]
  return(structure(list(
    data = mcaobj$data,
    memobj = mcaobj$memobj,
    Upyxcb = mcaobj$Upyxcb,
    test = list(permute = FALSE,
                significant = signif,
                test.table = ttable,
                details = NULL)),
    class = "mca"))
}
#
permute.mca <- function(mcaobj, permute = NA, alpha = 0.05, max.step = Inf) {
  permut <- function(p,ryx,us,t2,details,permute,step,ddfr){
     .Call("permut",p,ryx,us,
          t2,details,permute,
          step,ddfr, PACKAGE="codep")
  }
  if(!is.finite(max.step[1])) max.step <- ncol(mcaobj$memobj$U)
  if(is.na(permute[1])) permute <- minpermute(alpha,ncol(mcaobj$memobj$U),10,3)
  us <- matrix(NA, nrow(mcaobj$memobj$U), 0)
  uspyx <- matrix(NA, 0, 2)
  yxc <- cbind(yc=mcaobj$data[,1]-mean(mcaobj$data[,1]), xc=mcaobj$data[,2]-mean(mcaobj$data[,2]))
  ord <- order(abs(mcaobj$Upyxcb[,3]), decreasing = TRUE)
  ttable <- matrix(NA, 0, 4)
  colnames(ttable) <- c("t2", "ddf", "Testwise p", "Familywise p")
  details <- matrix(NA, 0, 3)
  colnames(details) <- c("t2* <= -|t2|", "-|t2| < t2* < |t2|", "t2* >= |t2|")
  step <- 1
  while(step != 0) {
    us <- cbind(us, mcaobj$memobj$U[,ord[step]])
    uspyx <- rbind(uspyx, mcaobj$Upyxcb[ord[step],1:2])
    ddfr <- nrow(mcaobj$data) - step - 1
    ryx <- yxc - (us %*% uspyx)
    t2 <- ddfr * (uspyx[step,1] * uspyx[step,2]) / sqrt((t(ryx[,1]) %*% ryx[,1]) * (t(ryx[,2]) %*% ryx[,2]))
    ttable <- rbind(ttable, c(t2, ddfr, NA, NA))
    details <- rbind(details,c(0,0,1))
    p <- matrix(NA, nrow(mcaobj$data), 4)
    details <- permut(p,ryx,us,t2,details,permute,step,ddfr)
    colnames(details) <- c("t2* <= -|t2|", "-|t2| < t2* < |t2|", "t2* >= |t2|")
    ttable[step,3] <- (details[step,1] + details[step,3]) / sum(details[step,])
    ttable[step,4] <- 1 - (1 - ttable[step,3])^(ncol(mcaobj$memobj$U) - step + 1)
    if (ttable[step,4] > alpha || step >= max.step) {
      rownames(ttable) <- colnames(mcaobj$memobj$U)[ord[1:step]]
      rownames(details) <- colnames(mcaobj$memobj$U)[ord[1:step]]
      step <- 0
    }
    else step <- step + 1
  }
  signif <- ord[which(ttable[,4] <= alpha)]
  return(structure(list(
    data = mcaobj$data,
    memobj = mcaobj$memobj,
    Upyxcb = mcaobj$Upyxcb,
    test = list(permute = permute,
                significant = signif,
                test.table = ttable,
                details = details)),
    class = "mca"))
}
#
print.mca <- function(x, ...) {
  cat("\nMulti-scale Codependence Analysis\n---------------------------------\n\n")
  cat("Coefficients:\n")
  print(cbind(round(x$Upyxcb[,3:4],5), Lambda=round(x$memobj$lambda,5)))
  cat("\n")
  return(invisible(NULL))
}
#
summary.mca <- function(object, ...) {
  if(is.null(object$test)) {
    cat("\nNo testing informations available\n\n")
  } else {
    cat("\nTest table:\n")
    print(object$test$test.table)
    cat("\n")
  }
  return(invisible(NULL))
}
#
plot.mca <- function(x, ...) {
  cc <- rep(grey(0.5),nrow(x$Upyxcb))
  if(!is.null(x$test)) cc[x$test$significant] <- grey(0)
  barplot(x$Upyxcb[,3],names.arg=rownames(x$Upyxcb),ylab="C",
          ylim=c(-1,1)*max(abs(x$Upyxcb[,3])), las=2, space=0, col = cc)
  return(invisible(NULL))
}
#
fitted.mca <- function(object, which=NA, components=FALSE, ...) {
  if(!is.null(object$test)) which <- object$test$significant
  else if(is.na(which[1])) stop("No testing informations available: user must identify relevant coefficients.")
  fit <- matrix(0,nrow(object$data),1)
  if(components) {
    cpns <- matrix(NA, nrow(object$data), length(which))
    rownames(cpns) <- rownames(object$memobj$U)
  }
  if(length(which) != 0) {
    by <- object$Upyxcb[which,4] * object$Upyxcb[which,2]
    fit <- object$memobj$U[,which] %*% cbind(by) + mean(object$data[,1])
    if (components) {
      for (i in 1:length(which)) cpns[,i] <- object$memobj$U[,which[i]] * by[i]
      colnames(cpns) <- paste("Component", which)
    }
  }
  colnames(fit) <- "fitted" ; rownames(fit) <- rownames(object$memobj$U)
  if(components) {
    return(list(fitted=fit, components=cpns))
  } else return(fit)
}
#
residuals.mca <- function(object, which=NA, ...) {
  if(!is.null(object$test)) which <- object$test$significant
  else if(is.na(which[1])) stop("No testing informations available: user must identify relevant coefficients.")
  res <- object$data[,1]
  if(length(which) != 0) {
    by <- object$Upyxcb[which,4] * object$Upyxcb[which,2]
    res <- res - object$memobj$U[,which] %*% cbind(by) - mean(object$data[,1])
  }
  colnames(res) <- "residuals" ; rownames(res) <- rownames(object$memobj$U)
  return(res)
}
#
predict.mca <- function(object, which=NA, newdata=NA, components=FALSE, ...) {
  if(is.na(newdata[1])) return (fitted.mca(object, which=which))
  if(!is.null(object$test)) which <- object$test$significant
  else if(is.na(which[1])) stop("No testing informations available: user must identify relevant coefficients.")
  if(is.matrix(newdata)) {
    newdata <- newdata[,1] ; warning("Only the first row of the matrix provided as 'newdata' is used.")
  }
  if(length(newdata) != nrow(object$memobj$U)) {
    stop("Number of observations in 'newdata' does not match the number of lines in U.")
  }
  mnew <- mean(newdata) ; newc <- cbind(newc=newdata-mnew)
  Upnew <- t(object$memobj$U) %*% newc
  pred <- matrix(0,nrow(object$data),1)
  if(components) {
    cpns <- matrix(NA, nrow(object$data), length(which))
    rownames(cpns) <- rownames(object$memobj$U)
  }
  if(length(which) != 0) {
    by <- object$Upyxcb[which,4] * Upnew[which]
    pred <- object$memobj$U[,which] %*% cbind(by) + mnew * mean(object$data[,1]) / mean(object$data[,2])
    if(components) {
      for (i in 1:length(which)) cpns[,i] <- object$memobj$U[,which[i]] * by[i]
      colnames(cpns) <- paste("Component", which)
      }
  }
  colnames(pred) <- "predicted" ; rownames(pred) <- rownames(object$memobj$U)
  if(components) {
    return(list(predicted=pred, components=cpns))
  } else return(pred)
}
#
