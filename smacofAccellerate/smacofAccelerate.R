library(MASS)
library(microbenchmark)
library(numDeriv)

source("smacofUtils.R")

smacofAccelerate <- function(delta,
                             ndim = 2,
                             wgth = 1 - diag(nrow(delta)),
                             xold = smacofTorgerson(delta, ndim),
                             opt = 1,
                             halt = 0,
                             wd = 4,
                             dg = 15,
                             itmax = 1000,
                             epsx = 1e-10,
                             epsf = 1e-15,
                             verbose = 1) {
  vmat <- -wgth
  diag(vmat) <- -rowSums(vmat)
  vinv <- ginv(vmat)
  nobj <- nrow(xold)
  xold <- xold %*% qr.Q(qr(t(xold[1:ndim, ])))
  bs <- smacofMakeBasis(nobj, ndim, wgth)
  cold <- Inf
  itel <- 1
  repeat {
    xold <- apply(xold, 2, function(x)
      x - mean(x))
    dold <- as.matrix(dist(xold))
    sold <- sum(wgth * (delta - dold) ^ 2)
    bold <- -wgth * delta / (dold + diag(nobj))
    diag(bold) <- -rowSums(bold)
    xbar <- vinv %*% bold %*% xold
    if (opt == 1) {
      h <- smacofOptionOne(xold, xbar, delta, wgth, vmat)
    }
    if (opt == 2) {
      h <- smacofOptionTwo(xold, xbar, delta, wgth, vmat)
    }
    if (opt == 3) {
      h <- smacofOptionThree(xold, xbar, delta, wgth, vmat)
    }
    if (opt == 4) {
      h <- smacofOptionFour(xold, xbar, delta, wgth, vmat, bs)
    }
    if (opt == 5) {
      h <- smacofOptionFive(xold, xbar, delta, wgth, vmat)
    }
    if (opt == 6) {
      h <- smacofOptionSix(xold, xbar, delta, wgth, vmat, vinv)
    }
    if (opt == 7) {
      h <- smacofOptionSeven(xold, xbar, delta, wgth, vmat)
    }
    if (opt == 8) {
      h <- smacofOptionEight(xold, xbar, delta, wgth, vmat, vinv)
    }
    labd <- sqrt(h$cnew / cold)
    if (verbose == 2) {
      smacofLinePrint(itel, sold, h$snew, h$cnew, labd, wd = wd, dg = dg)
    }
    if (halt == 1) {
      converge <- h$cnew < epsx
    } else {
      converge <- (sold - h$snew) < epsf
    }
    if ((itel == itmax) || converge) {
      break
    }
    itel <- itel + 1
    sold <- h$snew
    xold <- h$xnew
    cold <- h$cnew
  }
  if (verbose == 1) {
    smacofLinePrint(itel, sold, h$snew, h$cnew, labd, wd = wd, dg = dg)
  }
  if (opt == 5) {
    h$xnew <- (h$xnew + xold) / 2
    h$dnew <- as.matrix(dist(h$xnew))
    h$snew <- sum(wgth * (delta - h$dnew) ^ 2)
    smacofLinePrint(itel, sold, h$snew, h$cnew, labd, wd = wd, dg = dg)
  }
  if (opt == 6) {
    bold <- -wgth * delta / (h$dnew + diag(nobj))
    diag(bold) <- -rowSums(bold)
    h$xnew <- vinv %*% bold %*% h$xnew
    h$dnew <- as.matrix(dist(h$xnew))
    h$snew <- sum(wgth * (delta - h$dnew) ^ 2)
    smacofLinePrint(itel, sold, h$snew, h$cnew, labd, wd = wd, dg = dg)
  }
  return(
    list(
      x = h$xnew,
      s = h$snew,
      d = h$dnew,
      itel = itel,
      chng = h$cnew,
      labd = labd,
      wgth = wgth,
      delta = delta
    )
  )
}

smacofOptionOne <- function(xold, xbar, delta, wgth, vmat) {
  xnew <- xbar
  dnew <- as.matrix(dist(xnew))
  snew <- sum(wgth * (delta - dnew) ^ 2)
  cnew <- sum((xold - xnew) * (vmat %*% (xold - xnew)))
  return(list(
    xnew = xnew,
    dnew = dnew,
    snew = snew,
    cnew = cnew
  ))
}

smacofOptionTwo <- function(xold, xbar, delta, wgth, vmat) {
  ndim <- ncol(xold)
  xnew <- xbar %*% qr.Q(qr(t(xbar[1:ndim, ])))
  dnew <- as.matrix(dist(xnew))
  snew <- sum(wgth * (delta - dnew) ^ 2)
  cnew <- sum((xold - xnew) * (vmat %*% (xold - xnew)))
  return(list(
    xnew = xnew,
    dnew = dnew,
    snew = snew,
    cnew = cnew
  ))
}

smacofOptionThree <- function(xold, xbar, delta, wgth, vmat) {
  xnew <- xbar %*% svd(xbar)$v
  dnew <- as.matrix(dist(xnew))
  snew <- sum(wgth * (delta - dnew) ^ 2)
  cnew <- sum((xold - xnew) * (vmat %*% (xold - xnew)))
  return(list(
    xnew = xnew,
    dnew = dnew,
    snew = snew,
    cnew = cnew
  ))
}

smacofOptionFour <- function(xold, xbar, delta, wgth, vmat, bs) {
  ndim <- ncol(xold)
  nobj <- nrow(xold)
  xnew <- matrix(0, nobj, ndim)
  for (s in 1:ndim) {
    aux <- crossprod(bs[[s]], vmat %*% xbar[, s])
    xnew[, s] <- bs[[s]] %*% aux
  }
  dnew <- as.matrix(dist(xnew))
  snew <- sum(wgth * (delta - dnew) ^ 2)
  cnew <- sum((xold - xnew) * (vmat %*% (xold - xnew)))
  return(list(
    xnew = xnew,
    dnew = dnew,
    snew = snew,
    cnew = cnew
  ))
}

smacofOptionFive <- function(xold, xbar, delta, wgth, vmat) {
  xnew <- 2 * xbar - xold
  dnew <- as.matrix(dist(xnew))
  snew <- sum(wgth * (delta - dnew) ^ 2)
  cnew <- sum((xold - xnew) * (vmat %*% (xold - xnew)))
  return(list(
    xnew = xnew,
    dnew = dnew,
    snew = snew,
    cnew = cnew
  ))
}

smacofOptionSix <- function(xold, xbar, delta, wgth, vmat, vinv) {
  nobj <- nrow(xold)
  xaux <- 2 * xbar - xold
  daux <- as.matrix(dist(xaux))
  baux <- -wgth * delta / (daux + diag(nobj))
  diag(baux) <- -rowSums(baux)
  xbaz <- vinv %*% baux %*% xaux
  xnew <- 2 * xbaz - xaux
  dnew <- as.matrix(dist(xnew))
  snew <- sum(wgth * (delta - dnew) ^ 2)
  cnew <- sum((xold - xnew) * (vmat %*% (xold - xnew)))
  return(list(
    xnew = xnew,
    dnew = dnew,
    snew = snew,
    cnew = cnew
  ))
}

smacofOptionSeven <- function(xold, xbar, delta, wgth, vmat) {
  xaux <- 2 * xbar - xold
  daux <- as.matrix(dist(xaux))
  alpa <- sum(wgth * daux * delta) / sum(wgth * daux ^ 2)
  xnew <- alpa * xaux
  dnew <- alpa * daux
  snew <- sum(wgth * (delta - dnew) ^ 2)
  cnew <- sum((xold - xnew) * (vmat %*% (xold - xnew)))
  return(list(
    xnew = xnew,
    dnew = dnew,
    snew = snew,
    cnew = cnew
  ))
}

smacofOptionEight <- function(xold, xbar, delta, wgth, vmat, vinv) {
  nobj <- nrow(xold)
  xaux <- 2 * xbar - xold
  daux <- as.matrix(dist(xaux))
  baux <- -wgth * delta / (daux + diag(nobj))
  diag(baux) <- -rowSums(baux)
  xnew <- vinv %*% baux %*% xaux
  dnew <- as.matrix(dist(xnew))
  snew <- sum(wgth * (delta - dnew) ^ 2)
  cnew <- sum((xold - xnew) * (vmat %*% (xold - xnew)))
  return(list(
    xnew = xnew,
    dnew = dnew,
    snew = snew,
    cnew = cnew
  ))
}
