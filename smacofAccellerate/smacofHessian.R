numFunc <- function(x, nobj, ndim, wgth, delta, opt = 1) {
  xx <- matrix(x, nobj, ndim)
  dd <- as.matrix(dist(xx))
  vv <- -wgth
  diag(vv) <- -rowSums(vv)
  vinv <- solve(vv + (1 / nobj)) - (1 / nobj)
  bb <- -wgth * delta / (dd + diag(nobj))
  diag(bb) <- -rowSums(bb)
  xaux <- vinv %*% bb %*% xx
  if (opt == 1) {
    yy <- xaux
  }
  if (opt == 2) {
    lbd <- sqrt(sum(xaux[1, ] ^ 2))
    cs <- xaux[1, 2] / lbd
    sn <- xaux[1, 1] / lbd
    rot <- matrix(c(sn, cs, -cs, sn), 2, 2)
    yy <- xaux %*% rot
  }
  if (opt == 3) {
    yy <- xaux %*% svd(xaux)$v
  }
  if (opt == 4) {
    yy <- 2 * xaux - xx
  }
  return(as.vector(yy))
}

numHess <- function(x,
                    delta,
                    wgth = 1 - diag(nrow(x)),
                    opt = 1) {
  nobj <- nrow(x)
  ndim <- ncol(x)
  x <- as.vector(x)
  h <- jacobian(
    numFunc,
    x,
    nobj = nobj,
    ndim = ndim,
    wgth = wgth,
    delta = delta,
    opt = opt
  )
  return(h)
}

smacofRhoHessian <- function(x, delta, wgth) {
  n <- nrow(x)
  p <- ncol(x)
  np <- n * p
  dmat <- as.matrix(dist(x))
  fac1 <- wgth * delta / (dmat + diag(n))
  fac2 <- wgth * delta / ((dmat + diag(n)) ^ 3)
  bmat <- -fac1
  diag(bmat) <- -rowSums(bmat)
  hess <- matrix(0, np, np)
  for (s in 1:p) {
    ns <- (s - 1) * n + 1:n
    hess[ns, ns] <- bmat
    for (t in 1:p) {
      nt <- (t - 1) * n + 1:n
      ds <- outer(x[, s], x[, s], "-")
      dt <- outer(x[, t], x[, t], "-")
      aux <- -fac2 * ds * dt
      diag(aux) <- -rowSums(aux)
      hess[ns, nt] <- hess[ns, nt] - aux
    }
  }
  return(hess)
}

smacofJacobian <- function(x, delta, wgth) {
  n <- nrow(x)
  p <- ncol(x)
  np <- n * p
  vmat <- -wgth
  diag(vmat) <- -rowSums(vmat)
  vinv <- solve(vmat + (1 / n)) - (1 / n)
  jacob <- smacofRhoHessian(x, delta, wgth)
  for (s in 1:p) {
    ns <- (s - 1) * n + 1:n
    for (t in 1:p) {
      nt <- (t - 1) * n + 1:n
      jacob[ns, nt] <- vinv %*% jacob[ns, nt]
    }
  }
  return(jacob)
}
