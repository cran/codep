#
MCA <- function (Y, X, emobj) {
  if (!inherits(emobj, "eigenmap"))
    stop("Parameter 'emobj' must be a 'eigenmap' object!")
  if (!is.matrix(Y))
    Y <- matrix(Y, length(Y), 1L, dimnames = list(names(Y), "Y"))
  if (!is.matrix(X))
    X <- matrix(X, length(X), 1L, dimnames = list(names(X), "X"))
  if (nrow(Y) != nrow(X))
    stop("Number of observations in Y and X do not match!")
  if (nrow(emobj$U) != nrow(Y))
    stop("Number of observations in Y does not match the number of lines in U.")
  Dnms <- list(Y = colnames(Y), X = colnames(X))
  mssd <- list(mY = NA, mX = NA, ssdY = NA, ssdX = NA)
  mssd[[1L]] <- colMeans(Y)
  mssd[[2L]] <- colMeans(X)
  YXc <- list(Yc = Y - rep(mssd$mY, each = nrow(Y)),
              Xc = X - rep(mssd$mX, each = nrow(X)))
  mssd[[3L]] <- colSums(YXc$Yc^2)
  mssd[[4L]] <- colSums(YXc$Xc^2)
  UpYXcb <- list(UpY = matrix(NA, ncol(emobj$U), ncol(Y), dimnames = list(colnames(emobj$U), Dnms$Y)),
                 UpX = matrix(NA, ncol(emobj$U), ncol(X), dimnames = list(colnames(emobj$U), Dnms$X)),
                 C = array(NA, dim = c(ncol(emobj$U), ncol(Y), ncol(X)), dimnames = list(colnames(emobj$U), Dnms$Y, Dnms$X)),
                 B = array(NA, dim = c(ncol(emobj$U), ncol(Y), ncol(X)), dimnames = list(colnames(emobj$U), Dnms$Y, Dnms$X)),
                 CM = matrix(NA, ncol(emobj$U), ncol(X), dimnames = list(colnames(emobj$U), Dnms$X)))
  UpYXcb$UpY[] <- t(emobj$U) %*% YXc$Yc
  UpYXcb$UpX[] <- t(emobj$U) %*% YXc$Xc
  for (i in 1L:ncol(X)) {
    UpYXcb$C[,,i] <- UpYXcb$UpY/sqrt(rep(mssd[[3L]], each = ncol(emobj$U))) * UpYXcb$UpX[,i]/sqrt(mssd[[4L]][i])
    UpYXcb$B[,,i] <- UpYXcb$UpY/UpYXcb$UpX[,i]
    for (j in 1L:ncol(emobj$U)) {
      UpYXcb$CM[j,i] <- sqrt(sum((emobj$U[,j,drop=FALSE] %*% UpYXcb$UpY[j,,drop = FALSE])^2) / sum(YXc$Yc^2)) *
        abs(UpYXcb$UpX[j,i] / sqrt(mssd[[4L]][i]))
    }
  }
  return(structure(list(data = list(Y = Y, X = X), emobj = emobj, 
                        UpYXcb = UpYXcb, test = NULL), class = "cdp"))
}
#
test.cdp <- function(object, alpha = 0.05, max.step, response.tests = TRUE) {
  if (!inherits(object, "cdp"))
    stop("Parameter 'object' must be of class 'cdp'.")
  if (missing(max.step))
    max.step <- ncol(object$emobj$U)
  else
    max.step <- max.step[1L]
  #
  us <- matrix(NA, nrow(object$emobj$U), 1L)
  uspY <- matrix(NA, 1L, ncol(object$data$Y), dimnames = list(NULL, colnames(object$data$Y)))
  uspX <- matrix(NA, 1L, ncol(object$data$X), dimnames = list(NULL, colnames(object$data$X)))
  Yc <- object$data$Y - rep(colMeans(object$data$Y), each = nrow(object$data$Y))
  Xc <- object$data$X - rep(colMeans(object$data$X), each = nrow(object$data$X))
  ord <- order(apply(object$UpYXcb$CM, 1L, max), decreasing = TRUE)   # Order of eigenfunctions.
  bstX <- apply(object$UpYXcb$CM, 1L, which.max)                      # Which descriptor to test.
  ttable <- matrix(NA, 0L, 6L, dimnames = list(NULL, c("Variable","phi", "df1", "df2", "Testwise p", "Familywise p")))
  if (response.tests) 
    respts <- array(numeric(0), dim = c(ncol(object$data$Y), 4L, 0L),
                    dimnames = list(colnames(object$data$Y), c("tau", "df", "Testwise p", "Familywise p"), NULL))
  step <- 1L
  while (step != 0L) {
    us[] <- object$emobj$U[, ord[step]]
    uspY[] <- object$UpYXcb$UpY[ord[step],]
    uspX[] <- object$UpYXcb$UpX[ord[step],]
    df2 <- nrow(object$data$Y) - step - 1L
    Yhat <- us %*% uspY
    Xhat <- us %*% uspX
    Yc <- Yc - Yhat
    Xc <- Xc - Xhat
    phi_global <- df2^2 * sum(Yhat^2) * sum(Xhat[,bstX[ord[step]]]^2)/(sum(Yc^2) * sum(Xc[,bstX[ord[step]]]^2))
    ttable <- rbind(ttable, c(bstX[ord[step]], phi_global, ncol(object$data$Y), df2, NA, NA))
    ttable[step,5L] <- pphi(phi_global, ncol(object$data$Y), df2, lower.tail = FALSE)
    ttable[step,6L] <- 1 - (1 - ttable[step,5L])^((ncol(object$emobj$U) - step + 1) * ncol(object$data$X))
    if (response.tests) {
      tau_resp <- df2 * uspY[1L,] * uspX[bstX[ord[step]]] * (colSums(Yc^2) * sum(Xc[,bstX[ord[step]]]^2))^-0.5
      respts <- array(as.numeric(respts), dim = c(dim(respts)[1L], dim(respts)[2L], dim(respts)[3L] + 1L),
                      dimnames = c(dimnames(respts)[1L:2], list(NULL)))
      respts[,1L,step] <- tau_resp
      respts[,2L,step] <- df2
      respts[,3L,step] <- 2 * ptau(abs(tau_resp), df2, lower.tail = FALSE)
      respts[,4L,step] <- 1 - (1 - respts[,3L,step])^((ncol(object$emobj$U) - step + 1) * ncol(object$data$X))
    }
    if (ttable[step,6L] > alpha || step >= max.step) {
      rownames(ttable) <- colnames(object$emobj$U)[ord[1L:step]]
      if (response.tests)
        dimnames(respts)[[3L]] <- rownames(ttable)
      step <- 0
    }
    else step <- step + 1
  }
  signif <- list(U = ord[which(ttable[, 6L] <= alpha)])
  signif$X <- bstX[signif$U]
  return(structure(list(data = object$data, emobj = object$emobj, 
                        UpYXcb = object$UpYXcb, test = list(permute = FALSE, 
                                                            significant = signif, global = ttable, response = if (response.tests) respts else NULL, 
                                                            permutations = NULL)),
                   class = "cdp"))
}
#
permute.cdp <- function(object, permute, alpha = 0.05, max.step, response.tests = TRUE) {
  if (!inherits(object, "cdp"))
    stop("Parameter 'object' must be of class 'cdp'.")
  if (missing(max.step))
    max.step <- ncol(object$emobj$U)
  else
    max.step <- max.step[1L]
  if (missing(permute))
    permute <- minpermute(alpha, ncol(object$emobj$U) * ncol(object$data$X), 10L, 3L)
  else
    permute <- permute[1L]
  us <- matrix(NA, nrow(object$emobj$U), 1L)
  uspY <- matrix(NA, 1L, ncol(object$data$Y), dimnames = list(NULL, colnames(object$data$Y)))
  uspX <- matrix(NA, 1L, ncol(object$data$X), dimnames = list(NULL, colnames(object$data$X)))
  Yc <- object$data$Y - rep(colMeans(object$data$Y), each = nrow(object$data$Y))
  Xc <- object$data$X - rep(colMeans(object$data$X), each = nrow(object$data$X))
  ord <- order(apply(object$UpYXcb$CM, 1L, max), decreasing = TRUE)   # Order of eigenfunctions.
  bstX <- apply(object$UpYXcb$CM, 1L, which.max)                      # Which descriptor to test.
  ttable <- matrix(NA, 0L, 6L, dimnames = list(NULL, c("Variable","phi", "df1", "df2", "Testwise p", "Familywise p")))
  perm_global <- matrix(NA, 0L, 2L, dimnames = list(NULL, c("phi*<phi", "phi*>=phi")))
  if (response.tests) {
    respts <- array(numeric(0), dim = c(ncol(object$data$Y), 4L, 0L),
                    dimnames = list(colnames(object$data$Y), c("tau", "df", "Testwise p", "Familywise p"), NULL))
    perm_response <- array(numeric(0), dim = c(ncol(object$data$Y), 3L, 0L),
                           dimnames = list(colnames(object$data$Y), c("tau*<=-|tau|", "-|tau|<tau*<|tau|", "tau*>=|tau|"), NULL))
  }
  step <- 1L
  while (step != 0L) {
    us[] <- object$emobj$U[, ord[step]]
    uspY[] <- object$UpYXcb$UpY[ord[step],]
    uspX[] <- object$UpYXcb$UpX[ord[step],]
    df2 <- nrow(object$data$Y) - step - 1L
    Yhat <- us %*% uspY
    Xhat <- us %*% uspX
    Yc <- Yc - Yhat
    Xc <- Xc - Xhat
    phi_global0 <- sum(Yhat^2) * sum(Xhat[,bstX[ord[step]]]^2)/(sum(Yc^2) * sum(Xc[,bstX[ord[step]]]^2))
    ttable <- rbind(ttable, c(bstX[ord[step]], df2^2 * phi_global0, ncol(object$data$Y), df2, NA, NA))
    perm_global <- rbind(perm_global, c(0L, 1L))
    if (response.tests) {
      tau_resp0 <- uspY[1L,] * uspX[bstX[ord[step]]] * (colSums(Yc^2) * sum(Xc[, bstX[ord[step]]]^2))^-0.5
      respts <- array(as.numeric(respts), dim = c(dim(respts)[1L], dim(respts)[2L], dim(respts)[3L] + 1L),
                      dimnames = c(dimnames(respts)[1L:2], list(NULL)))
      respts[,1L,step] <- df2 * tau_resp0
      respts[,2L,step] <- df2
      perm_response <- array(as.numeric(perm_response),
                             dim = c(dim(perm_response)[1L], dim(perm_response)[2L], dim(perm_response)[3L] + 1L),
                             dimnames = c(dimnames(perm_response)[1L:2], list(NULL)))
      perm_response[,1L:2,step] <- 0
      perm_response[,3L,step] <- 1
      tmp <- .C("mcapermute",
                as.double(phi_global0),
                as.double(abs(tau_resp0)),
                as.double(Yc),
                as.integer(ncol(Yc)),
                as.double(Xc[,bstX[ord[step]]]),
                as.double(us),
                as.integer(nrow(us)),
                perm_global = as.integer(perm_global[step,]),
                perm_response = as.integer(perm_response[,,step]),
                as.integer(permute),
                as.integer(TRUE))
      perm_global[step,] <- tmp$perm_global
      perm_response[,,step] <- tmp$perm_response
    } else {
      perm_global[step,] <- .C("mcapermute",
                               as.double(phi_global0),
                               as.double(),
                               as.double(Yc),
                               as.integer(ncol(Yc)),
                               as.double(Xc[,bstX[ord[step]]]),
                               as.double(us),
                               as.integer(nrow(us)),
                               perm_global = as.integer(perm_global[step,]),
                               integer(),
                               as.integer(permute),
                               as.integer(FALSE))$perm_global
    }
    ttable[step,5L] <- perm_global[step,2L] / (permute + 1)
    ttable[step,6L] <- 1 - (1 - ttable[step, 5L])^((ncol(object$emobj$U) - step + 1) * ncol(object$data$X))
    if (response.tests) {
      respts[,3L,step] <- (perm_response[,1L,step] + perm_response[,3L,step]) / (permute + 1)
      respts[,4L,step] <- 1 - (1 - respts[,3L,step])^((ncol(object$emobj$U) - step + 1) * ncol(object$data$X))
    }
    if (ttable[step, 6L] > alpha || step >= max.step) {
      rownames(ttable) <- colnames(object$emobj$U)[ord[1L:step]]
      if (response.tests) {
        dimnames(respts)[[3L]] <- rownames(ttable)
        dimnames(perm_response)[[3L]] <- rownames(ttable)
      }
      step <- 0
    }
    else step <- step + 1
  }
  signif <- list(U = ord[which(ttable[,6L] <= alpha)])
  signif$X <- bstX[signif$U]
  return(structure(list(data = object$data, emobj = object$emobj, 
                        UpYXcb = object$UpYXcb, test = list(permute = permute, 
                                                            significant = signif, global = ttable, response = if (response.tests) respts else NULL, 
                                                            permutations = list(global = perm_global, response = if (response.tests) perm_response else NULL))),
                   class = "cdp"))
}
#
parPermute.cdp <- function (object, permute, alpha = 0.05, max.step, response.tests = TRUE, nnode, seeds, verbose = TRUE, ...) {
    if (!inherits(object, "cdp"))
        stop("Parameter 'object' must be of class 'cdp'.")
    if (missing(max.step))
        max.step <- ncol(object$emobj$U)
    else
        max.step <- max.step[1L]
    if (missing(permute))
        permute <- minpermute(alpha, ncol(object$emobj$U) * ncol(object$data$X), 10L, 3L)
    else
        permute <- permute[1L]
    if (missing(nnode))
        nnode <- detectCores()
    if(verbose)
        cat("Starting a cluster of", nnode, "nodes... ")
    cl <- makeCluster(nnode, ...)
    nnode <- length(cl)
    if(verbose)
        cat("done.\n")
    if(nnode < 2L)
        warning("Only a single worker could be recruited on that system.\nConsider using permute.mca() instead.")
    if(verbose)
        cat("Initializing workers:\n")
    if(missing(seeds))
        seeds <- as.integer(runif(nnode, -.Machine$integer.max, .Machine$integer.max))
    if(verbose)
        cat("Random seeds given to workers:", seeds, "\n")
    parSapply(cl = cl, X = seeds, FUN = function(x) set.seed(x))
    parSapply(cl = cl, X = 1L:nnode, function(x) require(codep))
    wpermute <- rep(permute%/%nnode, nnode)
    if(permute%%nnode)
        wpermute[1L:(permute%%nnode)] <- wpermute[1L:(permute%%nnode)]+1L
    C_mcapermute <- function(X, phi_global0,tau_ind0, rY, rx, us, perm_global, perm_response,ind)
        .C("mcapermute",
           as.double(phi_global0),
           as.double(abs(tau_ind0)),
           as.double(rY),
           as.integer(ncol(rY)),
           as.double(rx),
           as.double(us),
           as.integer(nrow(us)),
           as.integer(perm_global),
           as.integer(perm_response),
           as.integer(X),
           as.integer(ind))[8L:9L]
    us <- matrix(NA, nrow(object$emobj$U), 1L)
    uspY <- matrix(NA, 1L, ncol(object$data$Y), dimnames = list(NULL, colnames(object$data$Y)))
    uspX <- matrix(NA, 1L, ncol(object$data$X), dimnames = list(NULL, colnames(object$data$X)))
    Yc <- object$data$Y - rep(colMeans(object$data$Y), each = nrow(object$data$Y))
    Xc <- object$data$X - rep(colMeans(object$data$X), each = nrow(object$data$X))
    ord <- order(apply(object$UpYXcb$CM, 1L, max), decreasing = TRUE)
    bstX <- apply(object$UpYXcb$CM, 1L, which.max)
    ttable <- matrix(NA, 0L, 6L, dimnames = list(NULL, c("Variable", "phi", "df1", "df2", "Testwise p", "Familywise p")))
    perm_global <- matrix(NA, 0L, 2L, dimnames = list(NULL, c("phi*<phi", "phi*>=phi")))
    if (response.tests) {
        respts <- array(numeric(0), dim = c(ncol(object$data$Y), 4L, 0L),
                        dimnames = list(colnames(object$data$Y),
                                        c("tau", "df", "Testwise p", "Familywise p"), NULL))
        perm_response <- array(numeric(0), dim = c(ncol(object$data$Y), 3L, 0L),
                               dimnames = list(colnames(object$data$Y),
                                               c("tau*<=-|tau|", "-|tau|<tau*<|tau|", "tau*>=|tau|"),
                                               NULL))
    }
    if(verbose)
        cat("Performing permutation tests")
    step <- 1L
    while (step != 0L) {
        us[] <- object$emobj$U[,ord[step]]
        uspY[] <- object$UpYXcb$UpY[ord[step],]
        uspX[] <- object$UpYXcb$UpX[ord[step],]
        df2 <- nrow(object$data$Y) - step - 1L
        Yhat <- us %*% uspY
        Xhat <- us %*% uspX
        Yc <- Yc - Yhat
        Xc <- Xc - Xhat
        phi_global0 <- sum(Yhat^2) * sum(Xhat[,bstX[ord[step]]]^2)/(sum(Yc^2) * sum(Xc[,bstX[ord[step]]]^2))
        ttable <- rbind(ttable, c(bstX[ord[step]], df2^2 * phi_global0, ncol(object$data$Y), df2, NA, NA))
        perm_global <- rbind(perm_global, c(0L, 0L))
        if (response.tests) {
            tau_resp0 <- uspY[1L,] * uspX[bstX[ord[step]]] * (colSums(Yc^2) * sum(Xc[,bstX[ord[step]]]^2))^-0.5
            respts <- array(as.numeric(respts),
                            dim = c(dim(respts)[1L],dim(respts)[2L], dim(respts)[3L] + 1L),
                            dimnames = c(dimnames(respts)[1L:2L],list(NULL)))
            respts[, 1L, step] <- df2 * tau_resp0
            respts[, 2L, step] <- df2
            perm_response <- array(as.numeric(perm_response),
                                   dim = c(dim(perm_response)[1L], dim(perm_response)[2L],
                                           dim(perm_response)[3L] + 1L),
                                   dimnames = c(dimnames(perm_response)[1L:2L],list(NULL)))
            perm_response[,,step] <- 0
            wtmp <- parLapply(cl = cl, X = wpermute, fun = C_mcapermute, phi_global0 = phi_global0,
                              tau_ind0 = tau_resp0, rY = Yc, rx = Xc[,bstX[ord[step]]], us = us,
                              perm_global = perm_global[step,], perm_response = perm_response[,,step], ind=TRUE)
            tmp <- wtmp[[1L]]
            if(nnode > 1L)
                for(i in 2L:nnode) {
                    tmp[[1L]] <- tmp[[1L]] + wtmp[[i]][[1L]]
                    tmp[[2L]] <- tmp[[2L]] + wtmp[[i]][[2L]]
                }
            perm_global[step,] <- tmp[[1L]]
            perm_global[step,2L] <- perm_global[step,2L] + 1
            perm_response[,,step] <- tmp[[2L]]
            perm_response[,3L,step] <- perm_response[,3L,step] + 1
        }
        else {
            wtmp <- parLapply(cl = cl, X = wpermute, fun = C_mcapermute, phi_global0 = phi_global0,
                              tau_ind0 = tau_resp0, rY = Yc, rx = Xc[,bstX[ord[step]]], us = us,
                              perm_global = perm_global[step,], perm_response = integer(), ind = FALSE)
            tmp <- wtmp[[1L]][[1L]]
            if(nnode > 1L)
                for(i in 2L:nnode)
                    tmp <- tmp + wtmp[[i]][[1L]]
            perm_global[step,] <- tmp
            perm_global[step,2L] <- perm_global[step,2L] + 1
        }
        ttable[step,5L] <- perm_global[step,2L]/(permute + 1)
        ttable[step,6L] <- 1 - (1 - ttable[step,5L])^((ncol(object$emobj$U) - step + 1) * ncol(object$data$X))
        if (response.tests) {
            respts[,3L,step] <- (perm_response[,1L,step] + perm_response[,3L,step])/(permute + 1)
            respts[,4L,step] <- 1 - (1 - respts[,3L,step])^((ncol(object$emobj$U) - step + 1) * ncol(object$data$X))
        }
        if(verbose)
            cat(".")
        if (ttable[step,6L] > alpha || step >= max.step) {
            rownames(ttable) <- colnames(object$emobj$U)[ord[1L:step]]
            if (response.tests) {
                dimnames(respts)[[3L]] <- rownames(ttable)
                dimnames(perm_response)[[3L]] <- rownames(ttable)
            }
            step <- 0
        } else step <- step + 1
    }
    if(verbose)
        cat("done.\nStopping cluster.\n")
    parallel::stopCluster(cl)
    signif <- list(U = ord[which(ttable[, 6L] <= alpha)])
    signif$X <- bstX[signif$U]
    return(structure(list(data = object$data, emobj = object$emobj,
                          UpYXcb = object$UpYXcb,
                          test = list(permute = permute, significant = signif,
                                      global = ttable, response = if (response.tests) respts else NULL,
                                      permutations = list(global = perm_global,
                                                          response = if (response.tests) perm_response else NULL)),
                          seeds=seeds),
                     class = "cdp"))
}
#
print.cdp <- function (x, ...) {
  cat("\nMultiple Multi-scale Codependence Analysis\n---------------------------\n\n")
  cat(ncol(x$data$X), " explanatory variable", if (ncol(x$data$X) > 1L) "s","\n\n", sep = "")
  print(signif(cbind(x$emobj$lambda, x$UpYXcb$CM),4))
  if (!is.null(x$test)) {
    cat("--------------------\nGlobal testing information is available\n")
    if (!is.null(x$test$response)) 
      cat("Hypothesis test",if (ncol(x$data$X)>1L) "s"," also available for the response",
          if (ncol(x$data$X) > 1L) "s", "\n", sep = "")
    else cat("Individual test", if (ncol(x$data$X) > 1L) "s"," unavailable\n", sep = "")
  } else cat("\n")
  return(invisible(NULL))
}
#
plot.cdp <- function (x, col, col.signif = 2, main = "", ...) {
  if(missing(col))
    col <- grey(seq(1, 0, length.out = 256))
  mar <- par()$mar
  z <- log10(x$UpYXcb$CM + 1e-04)
  par(mar = c(mar[1L], mar[2L], mar[3L], 0.75), fig = c(0, 0.875 - 0.025 * (mar[4L] - 2.1), 0, 1))
  image(y = 1L:ncol(x$data$X), x = 1L:ncol(x$emobj$U), z = z,
        zlim = c(-4, 1e-04), col = col, axes = FALSE, xlab = "", ylab = "",
        main = main, ...)
  box(...)
  axis(1, at = 1L:ncol(x$emobj$U), labels = colnames(x$emobj$U), ...)
  axis(2, at = 1L:ncol(x$data$X), labels = colnames(x$data$X), ...)
  if (!is.null(x$test$signif)) 
    rect(xleft = x$test$signif$U - 0.5, xright = x$test$signif$U + 0.5,
         ybottom = x$test$signif$X - 0.5, ytop = x$test$signif$X + 0.5,
         border = col.signif, density = NULL, ...)
  par(mar = c(mar[1L], 0.75, mar[3L], mar[4L]), fig = c(0.875 - 0.025 * (mar[4L] - 2.1), 1, 0, 1), new = TRUE)
  image(z = matrix(seq(-4, 1e-04, length.out = 256), 1L, 256L), x = 0,
        y = seq(-4, 1e-04, length.out = 256), col = col, axes = FALSE,
        xlab = "", ylab = "", main = "", ...)
  box(...)
  axis(4, labels = 10^seq(-4, 0, 1), at = seq(-4, 0, 1), ...)
  par(mar = mar, fig = c(0, 1, 0, 1))
  return(invisible(NULL))
}
#
summary.cdp <- function(object, ...) {
  cat("\nMultiple Multi-scale Codependence Analysis\n---------------------------\n\n")
  cat(ncol(object$data$X), " explanatory variable",if (ncol(object$data$X)>1) "s","\n\n",sep="")
  if (is.null(object$test)) {
    cat("\nNo testing informations available\n\n")
  } else {
    cat("\nTest table:\n")
    tmp <- data.frame(object$test$global[,1L:4,drop=FALSE], character(nrow(object$test$global)),
                      stringsAsFactors = FALSE)
    colnames(tmp)[5L] <- "Pr(>|phi|)"
    tmp[,"Variable"] <- colnames(object$data$X)[tmp[,"Variable"]]
    tmp[object$test$global[,6L] < 2.2e-16,"Pr(>|phi|)"] <-"<2.2e-16"
    tmp[object$test$global[,6L] >= 2.2e-16,"Pr(>|phi|)"] <-
      signif(object$test$global[object$test$global[,6L] >= 2.2e-16,6L], 2L)
    print(tmp)
    cat("\n")
    if (is.null(object$test$response))
      cat("No individual response testing available\n")
    else {
      cat("Individual response tests:\n")
      tmp <- matrix(NA, nrow(object$test$global), 6L, dimnames = list(rownames(object$test$global), NULL))
      for (i in 1L:nrow(object$test$global))
        tmp[i,] <- hist(x = object$test$response[, 4L, i], plot = FALSE, breaks = c(0, 1e-04, 0.001, 0.01, 0.05, 0.1, 1))$counts
      tmp <- data.frame(tmp, Variable = colnames(object$data$X)[object$test$global[,1L]])
      colnames(tmp) <- c("p<=0.0001", "0.0001>p>=0.001", "0.001>p>=0.01", "0.01>p>=0.05", "0.05>p>=0.1", "p>0.1", "Variable")
      print(tmp)
      cat("\n")
    }
  }
  return(invisible(TRUE))
}
#
fitted.cdp <- function (object, selection, components = FALSE, ...) {
  if(missing(selection)) {
    if (!is.null(object$test))
      selection <- object$test$significant
    else
      stop("No testing informations available: user must identify the relevant coefficients.")
  } else {
    if (is.null(selection$U))
      stop("Parameter 'selection' must be a list with an element $U")
  }
  if (components) {
    cpns <- array(NA, dim = c(nrow(object$data$Y), ncol(object$data$Y), 
                              length(selection$U)))
    dimnames(cpns) <- list(rownames(object$data$Y), colnames(object$data$Y), 
                           colnames(object$emobj$U)[selection$U])
  }
  if (length(selection$U)) {
    by <- object$UpYXcb$UpY[selection$U,,drop=FALSE]
    fit <- object$emobj$U[,selection$U,drop=FALSE] %*% by + rep(colMeans(object$data$Y), each = nrow(object$data$Y))
    if (components) 
      for (i in 1L:length(selection$U))
        cpns[,,i] <- object$emobj$U[,selection$U[i],drop = FALSE] %*% by[i,,drop = FALSE]
  }
  if (components)
    return(list(fitted = fit, components = cpns))
  else
    return(fit)
}
#
residuals.cdp <- function (object, selection, ...) {
  if(missing(selection)) {
    if (!is.null(object$test))
      selection <- object$test$significant
    else
      stop("No testing informations available: user must identify the relevant coefficients.")
  } else {
    if (is.null(selection$U))
      stop("Parameter 'selection' must be a list with an element $U")
  }
  res <- object$data$Y - rep(colMeans(object$data$Y), each = nrow(object$data$Y))
  if (length(selection$U)) {
    by <- object$UpYXcb$UpY[selection$U,]
    res <- res - object$emobj$U[,selection$U,drop=FALSE] %*% by
  }
  return(res)
}
#
predict.cdp <- function (object, selection, newdata, components = FALSE, ...) {
  if (missing(newdata))
    return(fitted.cdp(object, selection = selection))
  if(missing(selection)) {
    if (!is.null(object$test))
      selection <- object$test$significant
    else
      stop("No testing informations available: user must identify relevant coefficients.")
  } else {
    if (is.null(selection$U)||is.null(selection$X))
      stop("Parameter 'selection' must be a list with elements $U and $X")
  }
  if(!is.null(newdata$X)) {
    if ((NROW(newdata$X) != nrow(object$data$X)) || (NCOL(newdata$X) != ncol(object$data$X)))
      stop("'newdata$X' (",NROW(newdata$X),"x",NCOL(newdata$X),") is not compatible with the original descriptors",
           nrow(object$data$X),"x",ncol(object$data$X))
    if (NROW(newdata$X) != nrow(object$emobj$U))
      stop("The number of observations in 'newdata$X' does not match the number of observations.")
    by <- matrix(NA, length(selection$U), ncol(object$data$Y),
                 dimnames = list(colnames(object$emobj$U)[selection$U], colnames(object$data$Y)))
    for (i in 1L:length(selection$X))
      by[i,] <- object$UpYXcb$B[selection$U[i],,selection$X[i]] * as.numeric(
        t(object$emobj$U[,selection$U[i],drop=FALSE]) %*% (newdata$X[,selection$X[i]] - mean(newdata$X[,selection$X[i]])))
  } else
    by <- object$UpYXcb$UpY[selection$U,,drop=FALSE]
  if (is.null(newdata$meanY))
    newdata$meanY <- colMeans(object$data$Y)
  else
    if (length(newdata$meanY) != ncol(object$data$Y))
      stop("The number of means in 'newdata$meanY' does not match the number of response variable.")
  if (is.null(newdata$target))
    newdata$target <- object$emobj$U
  else
    if (ncol(newdata$target) != ncol(object$emobj$U))
      stop("Incorrect number of eigenfunctions (columns) in 'newdata$meanY': ", ncol(newdata$target),
           ", while ",ncol(object$emobj$U)," is expected")
  if (components) {
    cpns <- array(NA, dim = c(nrow(newdata$target), ncol(object$data$Y), length(selection$U)))
    dimnames(cpns) <- list(rownames(newdata$target), colnames(object$data$Y), 
                           colnames(newdata$target)[selection$U])
  }
  pred <- newdata$target[,selection$U,drop=FALSE] %*% by + rep(newdata$meanY, each = nrow(newdata$target))
  if (components) {
    for (i in 1L:length(selection$U))
      cpns[,,i] <- newdata$target[,selection$U[i],drop=FALSE] %*% by[i,,drop=FALSE]
    return(list(predicted = pred, components = cpns))
  } else return(pred)
}
#
