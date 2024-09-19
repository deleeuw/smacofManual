
smacofUnweightWeighted <- function(delta,
                                   weights,
                                   ndim = 2,
                                   itmax = 1000,
                                   eps = 1e-15,
                                   verbose = TRUE) {
  nobj <- nrow(delta)
  weights <- weights / max(weights)
  xold <- smacofTorgerson(delta, ndim)
  dold <- as.matrix(dist(xold))
  sold <- sum(weights * (delta - dold) ^ 2)
  vmat <- -weights
  diag(vmat) <- -rowSums(vmat)
  itel <- 1
  repeat {
    bmat <- -weights * delta / (dold + diag(nobj))
    diag(bmat) <- -rowSums(bmat)
    xnew <- xold + ((bmat - vmat) %*% xold) / nobj
    dnew <- as.matrix(dist(xnew))
    snew <- sum(weights * (delta - dnew) ^ 2)
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
    itel <- itel + 1
  }
  return(list(
    x = xnew,
    s = snew,
    d = dnew,
    itel = itel
  ))
}

smacofWeighted <- function(delta,
                           weights,
                           ndim = 2,
                           itmax = 1000,
                           eps = 1e-15,
                           verbose = TRUE) {
  nobj <- nrow(delta)
  weights <- weights / max(weights)
  xold <- smacofTorgerson(delta, ndim)
  dold <- as.matrix(dist(xold))
  sold <- sum(weights * (delta - dold) ^ 2)
  vmat <- -weights
  diag(vmat) <- -rowSums(vmat)
  vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
  itel <- 1
  repeat {
    bmat <- -weights * delta / (dold + diag(nobj))
    diag(bmat) <- -rowSums(bmat)
    xnew <- vinv %*% (bmat %*% xold)
    dnew <- as.matrix(dist(xnew))
    snew <- sum(weights * (delta - dnew) ^ 2)
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
    itel <- itel + 1
  }
  return(list(
    x = xnew,
    s = snew,
    d = dnew,
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
