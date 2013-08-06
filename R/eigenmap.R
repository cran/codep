#
eigenmap <- function(x,opt.coord=NA,truncate=c(0,NA),wf,wpar,select=.Machine$double.eps^0.5) {
  #
  if(missing(wf)) stop("Parameter 'wf' must be one of 'binary', 'dbMEM', 'Drayf1', 'Drayf2', 'Drayf3', or 'sqrd'")
  if(!is.numeric(x)) stop("Parameter 'x' must be numeric!")
  wf <- match.arg(wf,c("binary","dbMEM","Drayf1","Drayf2","Drayf3","sqrd"))
  if(((wf=="Drayf2")||(wf=="Drayf3"))&&missing(wpar))
    stop("Parameter 'wpar' is mandatory with wf functions Drayf2 or Drayf3!")
  if(inherits(x,"dist")) {
    dd <- as.matrix(x)
    if(all(is.na(opt.coord))) {
      opt.coord <- as.matrix(x=1:nrow(dd))
      rownames(opt.coord) <- rownames(dd)
    } else {
      opt.coord <- as.matrix(opt.coord)
      if(nrow(opt.coord) != nrow(dd)) {
        stop(paste("You provided",nrow(opt.coord),"optional coordinates to reference", nrow(dd),"observations!"))
      } else rownames(opt.coord) <- rownames(dd)
    }
  } else {
    opt.coord <- as.matrix(x)
    x <- dist(x,method="euclidean") ; dd <- as.matrix(x)
    rownames(opt.coord) <- rownames(dd)    
  }
  if(wf[1] == "sqrd") {
    W <- -0.5*dd ; b <- NULL ; a <- NULL ; truncate <- NULL
  } else {
    if(any(is.na(truncate)) || length(truncate) != 2L) truncate <- c(0,max(hclust(x,method="single")$height))
    b <- matrix(0,nrow(dd),ncol(dd)) ; b[dd >= truncate[1L] & dd <= truncate[2]] <- 1
    if (wf == "binary") {
      W <- b ; a <- NULL
    } else {
      if (wf == "dbMEM") a <- 1-(dd/(4*truncate[2L]))^2
      if (wf == "Drayf1") a <- 1-dd/max(dd)
      if (wf == "Drayf2") a <- 1-(dd/max(dd))^wpar
      if (wf == "Drayf3") a <- 1/dd^wpar
      W <- b * a
    }
    diag(W) <- 0
  }
  n <- nrow(W) ; term <- diag(n) - (cbind(rep(1,n)) %*% rbind(rep(1,n)))/n ; O <- term %*% W %*% term
  eigO <- eigen(O)
  variables <- eigO$vectors[,abs(eigO$values) >= select]
  rownames(variables) <- rownames(dd) ; colnames(variables) <- paste("MEM",1L:ncol(variables),sep="")
  return(structure(list(
    coordinates=opt.coord,
    truncate=truncate,D=dd,weighting=wf,
    wpar=if((wf=="Drayf2")||(wf=="Drayf3")) wpar else NULL,
    lambda=eigO$values[abs(eigO$values) >= select],
    U=variables),class="eigenmap"))
}
#
print.eigenmap <- function(x, ...) {
  cat(paste("\nMoran's eigenvector map containing",length(x$lambda),"basis functions.\n"))
  cat(paste("Functions span",nrow(x$U),"observations.\n\n"))
  cat(paste("Eigenvalues:\n"))
  print.default(x$lambda)
  cat("\n")
  return(invisible(NULL))
}
#
plot.eigenmap <- function(x, ...) {
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
