smacofRobustPseudoHuber <- function(delta,
                                    weights = 1 - diag(nrow(delta)),
                                    ndim = 2,
                                    cons = 0,
                                    itmax = 1000,
                                    eps = 1e-15,
                                    verbose = TRUE) {
  nobj <- nrow(delta)
  wmax <- max(weights)
  xold <- smacofTorgerson(delta, ndim)
  dold <- as.matrix(dist(xold))
  rold <- sqrt((delta - dold) ^ 2 + cons)
  sold <- sum(weights * rold)
  wold <- weights / (rold + diag(nobj))
  itel <- 1
  repeat {
    vmat <- -wold
    diag(vmat) <- -rowSums(vmat)
    vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
    bmat <- -wold * delta / (dold + diag(nobj))
    diag(bmat) <- -rowSums(bmat)
    xnew <- vinv %*% (bmat %*% xold)
    dnew <- as.matrix(dist(xnew))
    rnew <- sqrt((delta - dnew) ^ 2 + cons)
    wnew <- weights / (rnew + diag(nobj))
    snew <- sum(weights * rnew)
    if (verbose) {
      cat(
        "itel ",
        formatC(itel, width = 4, format = "d"),
        "sold ",
        formatC(sold, digits = 10, format = "f"),
        "snew ",
        formatC(snew, digits = 10, format = "f"),
        "\n"
      )
    }
    if ((itel == itmax) || ((sold - snew) < eps)) {
      break
    }
    xold <- xnew
    dold <- dnew
    sold <- snew
    wold <- wnew
    rold <- rnew
    itel <- itel + 1
  }
  return(list(
    x = xnew,
    s = snew,
    d = dnew,
    r = rnew,
    itel = itel
  ))
}

smacofRobustHuber <- function(delta,
                              weights = 1 - diag(nrow(delta)),
                              ndim = 2,
                              cons = 0,
                              itmax = 1000,
                              eps = 1e-10,
                              verbose = TRUE) {
  nobj <- nrow(delta)
  wmax <- max(weights)
  xold <- smacofTorgerson(delta, ndim)
  dold <- as.matrix(dist(xold))
  fold <- abs(delta - dold)
  rold <- ifelse(fold < cons, (fold ^ 2) / 2, cons * fold - (cons ^ 2) / 2)
  sold <- sum(weights * rold)
  wold <- ifelse(fold < cons, weights, cons * weights / (fold + diag(nobj)))
  itel <- 1
  repeat {
    vmat <- -wold
    diag(vmat) <- -rowSums(vmat)
    vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
    bmat <- -wold * delta / (dold + diag(nobj))
    diag(bmat) <- -rowSums(bmat)
    xnew <- vinv %*% (bmat %*% xold)
    dnew <- as.matrix(dist(xnew))
    fnew <- abs(delta - dnew)
    rnew <- ifelse(fnew < cons, (fnew ^ 2) / 2, cons * fnew - (cons ^ 2) / 2)
    snew <- sum(weights * rnew)
    wnew <- ifelse(fnew < cons, weights, cons * weights / (fnew + diag(nobj)))
    if (verbose) {
      cat(
        "itel ",
        formatC(itel, width = 4, format = "d"),
        "sold ",
        formatC(sold, digits = 10, format = "f"),
        "snew ",
        formatC(snew, digits = 10, format = "f"),
        "\n"
      )
    }
    if ((itel == itmax) || ((sold - snew) < eps)) {
      break
    }
    xold <- xnew
    dold <- dnew
    sold <- snew
    wold <- wnew
    rold <- rnew
    itel <- itel + 1
  }
  return(list(
    x = xnew,
    s = snew,
    d = dnew,
    r = rnew,
    itel = itel
  ))
}

smacofRobustTukey <- function(delta,
                              weights = 1 - diag(nrow(delta)),
                              ndim = 2,
                              cons = 0,
                              itmax = 1000,
                              eps = 1e-10,
                              verbose = TRUE) {
  nobj <- nrow(delta)
  wmax <- max(weights)
  xold <- smacofTorgerson(delta, ndim)
  dold <- as.matrix(dist(xold))
  fold <- delta - dold
  rold <- ((cons ^ 2) / 6) * ifelse(abs(fold) < cons, (1 - (1 - (fold / cons) ^ 2) ^ 3), 1)
  sold <- sum(weights * rold)
  wold <- ifelse(abs(fold) < cons, weights * (1 - (fold / cons) ^ 2) ^ 2, 0) / 2
  itel <- 1
  repeat {
    vmat <- -wold
    diag(vmat) <- -rowSums(vmat)
    vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
    bmat <- -wold * delta / (dold + diag(nobj))
    diag(bmat) <- -rowSums(bmat)
    xnew <- vinv %*% (bmat %*% xold)
    dnew <- as.matrix(dist(xnew))
    fnew <- delta - dnew
    rnew <- ((cons ^ 2) / 6) * ifelse(abs(fnew) < cons, (1 - (1 - (fnew / cons) ^ 2) ^ 3), 1)
    snew <- sum(weights * rnew)
    wnew <- ifelse(abs(fnew) < cons, weights * (1 - (fnew / cons) ^ 2) ^ 2, 0) / 2
    if (verbose) {
      cat(
        "itel ",
        formatC(itel, width = 4, format = "d"),
        "sold ",
        formatC(sold, digits = 10, format = "f"),
        "snew ",
        formatC(snew, digits = 10, format = "f"),
        "\n"
      )
    }
    if ((itel == itmax) || ((sold - snew) < eps)) {
      break
    }
    xold <- xnew
    dold <- dnew
    sold <- snew
    wold <- wnew
    rold <- rnew
    itel <- itel + 1
  }
  return(list(
    x = xnew,
    s = snew,
    d = dnew,
    r = rnew,
    itel = itel
  ))
}

smacofRobustConvolution <- function(delta,
                              weights = 1 - diag(nrow(delta)),
                              ndim = 2,
                              cons = 0,
                              itmax = 1000,
                              eps = 1e-10,
                              verbose = TRUE) {
  nobj <- nrow(delta)
  wmax <- max(weights)
  xold <- smacofTorgerson(delta, ndim)
  dold <- as.matrix(dist(xold))
  fold <- delta - dold
  rold <- fold * (2 * pnorm(fold / cons) - 1) + 2 * cons * dnorm(fold / cons)
  sold <- sum(weights * rold)
  wold <- (pnorm(fold / cons) - 0.5) / (fold + diag(nobj))
  itel <- 1
  repeat {
    vmat <- -wold
    diag(vmat) <- -rowSums(vmat)
    vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
    bmat <- -wold * delta / (dold + diag(nobj))
    diag(bmat) <- -rowSums(bmat)
    xnew <- vinv %*% (bmat %*% xold)
    dnew <- as.matrix(dist(xnew))
    fnew <- delta - dnew
    rnew <- fnew * (2 * pnorm(fnew / cons) - 1) + 2 * cons * dnorm(fnew / cons)
    snew <- sum(weights * rnew)
    wnew <- (pnorm(fnew / cons) - 0.5) / (fnew + diag(nobj))
    if (verbose) {
      cat(
        "itel ",
        formatC(itel, width = 4, format = "d"),
        "sold ",
        formatC(sold, digits = 10, format = "f"),
        "snew ",
        formatC(snew, digits = 10, format = "f"),
        "\n"
      )
    }
    if ((itel == itmax) || ((sold - snew) < eps)) {
      break
    }
    xold <- xnew
    dold <- dnew
    sold <- snew
    wold <- wnew
    rold <- rnew
    itel <- itel + 1
  }
  return(list(
    x = xnew,
    s = snew,
    d = dnew,
    r = rnew,
    itel = itel
  ))
}

smacofTorgerson <- function(delta, ndim) {
  dd <- delta ^ 2
  rd <- apply(dd, 1, mean)
  md <- mean(dd)
  sd <- -.5 * (dd - outer(rd, rd, "+") + md)
  ed <- eigen(sd)
  return(ed$vectors[, 1:ndim] %*% diag(sqrt(ed$values[1:ndim])))
}
