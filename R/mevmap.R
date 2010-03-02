#
### Computation of Moran's eigenvector maps.
### Calculations outlined in Dray et al. (2006) Ecological Modelling 196: 483-493
### meigmap() returns a list object with a 'mem' class attribute.
### Dedicated printing and plotting functions are provided for class 'mem'
### Required package: stats
###
### Guillaume Guénard, January 2009
#
meigmap <- function(x,
                    opt.coordinates=NA,
                    truncation=c(0,NA),
                    weighting=c("1","f1","f2","f3"),
                    wpar=1,
                    select=1e-9) {
  #
  if(!is.numeric(x)) stop("Parameter 'x' must be numeric!")
  if(inherits(x,"dist")) {
    dd <- as.matrix(x)
    if(all(is.na(opt.coordinates))) {
      opt.coordinates <- as.matrix(x=1:nrow(dd))
      rownames(opt.coordinates) <- rownames(dd)
    } else {
      opt.coordinates <- as.matrix(opt.coordinates)
      if(nrow(opt.coordinates) != nrow(dd)) {
        stop(paste("You provided",nrow(opt.coordinates),"optional coordinates to reference", nrow(dd),"observations!"))
      } else rownames(opt.coordinates) <- rownames(dd)
    }
  } else {
    opt.coordinates <- as.matrix(x)
    x <- dist(x,method="euclidean") ; dd <- as.matrix(x)
    rownames(opt.coordinates) <- rownames(dd)    
  }
  if (any(is.na(truncation)) || length(truncation) != 2) truncation <- c(0,max(hclust(x,method="single")$height))
  b <- matrix(0,nrow(dd),ncol(dd)) ; b[dd > truncation[1] & dd <= truncation[2]] <- 1
  if (!any(weighting[1] == c("1","f1","f2","f3"))) stop("Parameter 'weighting' must be either 1, f1, f2, or f3!")
  if (weighting[1] == 1) a <- 1
  if (weighting[1] == "f1") a <- 1 - dd/max(dd)
  if (weighting[1] == "f2") a <- 1 - (dd/max(dd))^wpar
  if (weighting[1] == "f3") {a <- 1/dd^wpar ; diag(a) <- 0}
  W <- b * a
  n <- nrow(W) ; term <- diag(n) - (cbind(rep(1,n)) %*% rbind(rep(1,n)))/n ; O <- term %*% W %*% term
  eigO <- eigen(O)
  variables <- eigO$vectors[,abs(eigO$values) >= select]
  rownames(variables) <- rownames(dd) ; colnames(variables) <- paste("MEM",1:ncol(variables),sep="")
  return(structure(list(
    coordinates=opt.coordinates,
    lambda=eigO$values[abs(eigO$values) >= select],
    U=variables),class="mem"))
}
#
mem <- function(coordinates,lambda,U) {
  # Verify that the parameters are input correctly.
  if(!is.numeric(coordinates)) stop("Parameter 'coordinate' must be of type 'numeric'!")
  if(!is.numeric(lambda)) stop("Parameter 'lambda' must be of type 'numeric'!")
  if(!is.numeric(U)) stop("Parameter 'U' must be of type 'numeric'!")
  if(!is.matrix(coordinates)) coordinates <- as.matrix(coordinates)
  if(is.matrix(lambda)) stop("Provide the parameter 'lambda' (eigenvalues) as a numeric vector, not a matrix!")
  if(!is.matrix(U)) stop("Parameter 'U' must be a matrix!")
  if(nrow(coordinates) != nrow(U)) stop("You provided",nrow(coordinates),"coordinates for",nrow(U),"observations!")
  if(length(lambda) != ncol(U)) stop("You provided",nrow(coordinates),"eigenvalues for",ncol(U),"eigenvectors!")
  return(structure(list(coordinates=coordinates, lambda=lambda, U=U), class="mem"))
}
#
print.mem <- function(x, ...) {
  cat(paste("\nMoran's eigenvector map containing",length(x$lambda),"basis functions.\n"))
  cat(paste("Functions span",nrow(x$U),"observations.\n\n"))
  cat(paste("Eigenvalues:\n"))
  print.default(x$lambda)
  cat("\n")
  return(invisible(NULL))
}
#
plot.mem <- function(x, ...) {
  if (ncol(x$coordinates) > 2) {
    warning(paste(ncol(x$coordinates),"dimensions were provided but only the first 2 were used for plotting."))
  }
  cat("Left-click on the graphical display to see further variables or right-click (Mac: esc) to terminate plotting.\n")
  for (i in 1:length(x$lambda)) {
    layout(matrix(c(1,2,2),1,3))
    plot(y=x$lambda,x=1:length(x$lambda),ylab=expression(lambda),xlab="Order",type="b")
    title("Eigenvalues diagram")
    points(y=x$lambda[i],x=i,pch=21,bg="black")
    if (ncol(x$coordinates) == 1) plot(y=x$U[,i],x=x$coordinates[,1],ylab="value",xlab="Location",type="l")
    if (ncol(x$coordinates) > 1) {
      plot(y=x$coordinates[,2],x=x$coordinates[,1],
           asp=1,ylab="Location (y)",xlab="Location (x)",type="p",pch=3)
      gcol <- grey((sign(x$U[,i])+1)/2)
      gsize <- 3*abs(x$U[,i])/max(abs(x$U[,i]))
      points(y=x$coordinates[,2],x=x$coordinates[,1], pch = 21, bg=gcol, cex = gsize)
    }
    title(paste("Variable display:",colnames(x$U)[i]))
    ttt <- locator(1)
    if (is.null(ttt)) break
  }
  return(invisible(NULL))
}
#
