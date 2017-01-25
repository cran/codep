#
cthreshold <- function(alpha,nbtest) return(1 - (1 - alpha)^(nbtest^-1))
#
minpermute <- function(alpha,nbtest,margin=1,ru=3) {
  return(round(floor(margin * (1 - (1 - alpha)^(nbtest^-1))^-1), -ru) + (10^ru) - 1)
}
#
dphi <- function(x, nu1, nu2, tol = .Machine$double.eps ^ 0.5) {
  res <- numeric(length(x))
  nu1 <- rep(nu1, length.out = length(x))
  nu2 <- rep(nu2, length.out = length(x))
  # Integrand function to be integrated over positive real numbers.
  f <- function(z, x, nu1, nu2)
    df(z, nu1, nu1 * nu2) * df(x / z, 1, nu2) / z   # abs(z) always == z for real positive.
  # The integrand is strictly positive.
  for(i in 1L:length(x))
    res[i] <- integrate(f, lower = 0, upper = Inf, x = x[i], nu1 = nu1[i], nu2 = nu2[i], rel.tol = tol)$value
  res
}
#
pphi <- function(q, nu1, nu2, lower.tail = TRUE, tol = .Machine$double.eps^0.5) {
  res <- numeric(length(q))
  nu1 <- rep(nu1, length.out = length(q))
  nu2 <- rep(nu2, length.out = length(q))
  if(lower.tail)
    for(i in 1L:length(q))
      res[i] <- 1 - integrate(dphi, lower = q[i], upper = Inf, nu1 = nu1[i], nu2 = nu2[i], tol = tol, rel.tol = tol)$value
  else
    for(i in 1L:length(q))
      res[i] <- integrate(dphi, lower = q[i], upper = Inf, nu1 = nu1[i], nu2 = nu2[i], tol = tol, rel.tol = tol)$value
  res
}
#
dtau <- function(x, nu, tol = .Machine$double.eps ^ 0.5) {
  res <- numeric(length(x))
  nu <- rep(nu, length.out = length(x))
  # Integrand function to be integrated over positive real numbers.
  f <- function(z, x, nu)
    dt(z, nu) * dt(x / z, nu) / z   # abs(z) always == z for real positive.
  # The integrand being symmetric, only the positive is integrated and the result multiplied by two.
  for(i in 1L:length(x))
    res[i] <- 2 * integrate(f, lower = 0, upper = Inf, x = x[i], nu = nu[i], rel.tol = tol)$value
  res
}
#
ptau <- function(q, nu, lower.tail = TRUE, tol = .Machine$double.eps^0.5) {
  res <- rep(0.5, length(q))
  nu <- rep(nu, length.out = length(q))
  # Only the positive 
  signq <- sign(q)
  q <- abs(q)
  # Code avoid calculate PDF(0) because the function has a singularity at that point.
  if(lower.tail) {
    for(i in which(!!signq)) {
      if(signq[i]==1)
        res[i] <- 1 - integrate(dtau, lower = q[i], upper = Inf, nu = nu[i], tol = tol, rel.tol = tol)$value
      else
        res[i] <- integrate(dtau, lower = q[i], upper = Inf, nu = nu[i], tol = tol, rel.tol = tol)$value
    }
  } else {
    for(i in which(!!signq)) {
      if(signq[i]==1)
        res[i] <- integrate(dtau, lower = q[i], upper = Inf, nu = nu[i], tol = tol, rel.tol = tol)$value
      else
        res[i] <- 1 - integrate(dtau, lower = q[i], upper = Inf, nu = nu[i], tol = tol, rel.tol = tol)$value
    }
  }
  res
}
#
