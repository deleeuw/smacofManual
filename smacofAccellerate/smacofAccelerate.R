library(MASS)
library(microbenchmark)
library(numDeriv)

smacofAccelerate <- function(delta,
                             ndim = 2,
                             wgth = 1 - diag(nrow(delta)),
                             xold = smacofTorgerson(delta, ndim),
                             opt = 1,
                             halt = 1,
                             itmax = 1000,
                             epsx = 1e-10,
                             epsf = 1e-15,
                             verbose = 1) {
  v <- -wgth
  diag(v) <- -rowSums(v)
  vinv <- ginv(v)
  n <- nrow(xold)
  cold <- Inf
  itel <- 1
  last <- FALSE
  repeat {
    xold <- apply(xold, 2, function(x)
      x - mean(x))
    dold <- as.matrix(dist(xold))
    sold <- sum(wgth * (delta - dold) ^ 2)
    bold <- -wgth * delta / (dold + diag(n))
    diag(bold) <- -rowSums(bold)
    xbar <- vinv %*% bold %*% xold
    if (opt == 1) {
      xnew <- xbar
      dnew <- as.matrix(dist(xnew))
      snew <- sum(wgth * (delta - dnew) ^ 2)
      cnew <- sum((xold - xnew) * (v %*% (xold - xnew)))
    }
    if (opt == 2) {
      lbd <- sqrt(sum(xbar[1, ] ^ 2))
      cs <- xbar[1, 2] / lbd
      sn <- xbar[1, 1] / lbd
      rot <- matrix(c(sn, cs, -cs, sn), 2, 2)
      xnew <- xbar %*% rot
      dnew <- as.matrix(dist(xnew))
      snew <- sum(wgth * (delta - dnew) ^ 2)
      cnew <- sum((xold - xnew) * (v %*% (xold - xnew)))
    }
    if (opt == 3) {
      xnew <- xbar %*% svd(xbar)$v
      dnew <- as.matrix(dist(xnew))
      snew <- sum(wgth * (delta - dnew) ^ 2)
      cnew <- sum((xold - xnew) * (v %*% (xold - xnew)))
    }
    if (opt == 4) {
      xnew <- 2 * xbar - xold
      dnew <- as.matrix(dist(xnew))
      snew <- sum(wgth * (delta - dnew) ^ 2)
      cnew <- sum((xold - xnew) * (v %*% (xold - xnew)))
    }
    if (opt == 5) {
      xaux <- 2 * xbar - xold
      daux <- as.matrix(dist(xaux))
      baux <- -wgth * delta / (daux + diag(n))
      diag(baux) <- -rowSums(baux)
      xbaz <- vinv %*% baux %*% xaux
      xnew <- 2 * xbaz - xaux
      dnew <- as.matrix(dist(xnew))
      snew <- sum(wgth * (delta - dnew) ^ 2)
      cnew <- sum((xold - xnew) * (v %*% (xold - xnew)))
    }
    if (opt == 6) {
      xaux <- 2 * xbar - xold
      daux <- as.matrix(dist(xaux))
      alpa <- sum(wgth * daux * delta) / sum(wgth * daux ^ 2)
      xnew <- alpa * xaux
      dnew <- alpa * daux
      snew <- sum(wgth * (delta - dnew) ^ 2)
      cnew <- sum((xold - xnew) * (v %*% (xold - xnew)))
    }
    if (opt == 7) {
      xaux <- 2 * xbar - xold
      daux <- as.matrix(dist(xaux))
      baux <- -wgth * delta / (daux + diag(n))
      diag(baux) <- -rowSums(baux)
      xnew <- vinv %*% baux %*% xaux
      dnew <- as.matrix(dist(xnew))
      snew <- sum(wgth * (delta - dnew) ^ 2)
      cnew <- sum((xold - xnew) * (v %*% (xold - xnew)))
    }
    labd <- sqrt(cnew / cold)
    if (verbose == 2) {
      cat(
        "itel",
        formatC(itel, digits = 2, format = "d"),
        "sold",
        formatC(sold, digits = 10, format = "f"),
        "snew",
        formatC(snew, digits = 10, format = "f"),
        "chng",
        formatC(cnew, digits =  10, format = "f"),
        "labd",
        formatC(labd, digits =  10, format = "f"),
        "\n"
      )
    }
    if (halt == 1) {
      converge <- cnew < epsx
    } else {
      converge <- (sold - snew) < epsf
    }
    if ((itel == itmax) || converge) {
      break
    }
    itel <- itel + 1
    sold <- snew
    xold <- xnew
    cold <- cnew
  }
  if (verbose == 1) {
    cat(
      "itel",
      formatC(itel, digits = 2, format = "d"),
      "sold",
      formatC(sold, digits = 10, format = "f"),
      "snew",
      formatC(snew, digits = 10, format = "f"),
      "chng",
      formatC(cnew, digits =  10, format = "f"),
      "labd",
      formatC(labd, digits =  10, format = "f"),
      "\n"
    )
  }
  if (opt == 4) {
    xaux <- (xnew + xold) / 2
    xnew <- (xnew + xold) / 2
    dnew <- as.matrix(dist(xnew))
    snew <- sum(wgth * (delta - dnew) ^ 2)
  }
  if (opt == 5) {
    bold <- -wgth * delta / (dnew + diag(n))
    diag(bold) <- -rowSums(bold)
    xnew <- vinv %*% bold %*% xnew
    dnew <- as.matrix(dist(xnew))
    snew <- sum(wgth * (delta - dnew) ^ 2)
  }
  return(list(
    x = xnew,
    s = snew,
    d = dnew,
    itel = itel,
    chng = cnew,
    labd = labd,
    wgth = wgth,
    delta = delta
  ))
}

smacofCompare <- function(delta, ndim = 2) {
  n <- nrow(delta)
  wgth <- 1 - diag(n)
  xold <- smacofTorgerson(delta, ndim)
  return(
    microbenchmark(
      smacofAccelerate(
        delta,
        ndim = 2,
        opt = 1,
        halt = 2,
        verbose = FALSE
      ),
      smacofAccelerate(
        delta,
        ndim = 2,
        opt = 2,
        halt = 2,
        verbose = FALSE
      ),
      smacofAccelerate(
        delta,
        ndim = 2,
        opt = 3,
        halt = 2,
        verbose = FALSE
      ),
      smacofAccelerate(
        delta,
        ndim = 2,
        opt = 4,
        halt = 2,
        verbose = FALSE
      ),
      smacofAccelerate(
        delta,
        ndim = 2,
        opt = 5,
        halt = 2,
        verbose = FALSE
      ),
      smacofAccelerate(
        delta,
        ndim = 2,
        opt = 6,
        halt = 2,
        verbose = FALSE
      ),
      smacofAccelerate(
        delta,
        ndim = 2,
        opt = 7,
        halt = 2,
        verbose = FALSE
      )
    )
  )
}

smacofTorgerson <- function(delta, ndim) {
  n <- nrow(delta)
  dd <- delta ^ 2
  rd <- rowSums(dd) / n
  sd <- mean(dd)
  cc <- -.5 * (dd - outer(rd, rd, "+") + sd)
  ee <- eigen(cc)
  x <- ee$vectors[, 1:ndim] %*% diag(sqrt(ee$values[1:ndim]))
  return(x)
}

numFunc <- function(x, nobj, ndim, wgth, delta, opt = 1) {
  xx <- matrix(x, nobj, ndim)
  dd <- as.matrix(dist(xx))
  vv <- -wgth
  diag(vv) <- -rowSums(vv)
  vinv <- ginv(vv)
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

numHess <- function(x, delta, wgth = 1 - diag(nrow(x)), opt = 1) {
  nobj <- nrow(x)
  ndim <- ncol(x)
  x <- as.vector(x)
  h <- jacobian(numFunc, x, nobj = nobj, ndim = ndim, wgth = wgth, delta = delta, opt = opt )
  return(h)
}
