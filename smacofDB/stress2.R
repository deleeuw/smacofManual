
stress2 <-
  function(delta,
           wmat = 1 - diag(nrow(delta)),
           ndim = 2,
           itmax = 1000,
           eps = 1e-10,
           verbose = TRUE) {
    itel <- 1
    n <- nrow(delta)
    wmat <- wmat / sum(wmat)
    vmat <- -wmat
    diag(vmat) <- -rowSums(vmat)
    xold <- torgerson(delta, ndim)
    dold <- as.matrix(dist(xold))
    enum <- sum(wmat * delta * dold)
    eden <- sum(wmat * dold ^ 2)
    lbda <- enum / eden
    dold <- lbda * dold
    xold <- lbda * xold
    aold <- sum(wmat * dold)
    sold <- sum(wmat * (delta - dold) ^ 2) / sum(wmat * (dold - aold) ^ 2)
    repeat {
      mmat <- -aold * wmat / (dold + diag(n))
      diag(mmat) <- -rowSums(mmat)
      bmat <- -wmat * delta / (dold + diag(n))
      diag(bmat) <- -rowSums(bmat)
      umat <- ((1 - sold) * vmat) + (sold * mmat)
      uinv <- solve(umat + 1/n) - 1/n
      xnew <- uinv %*% bmat %*% xold
      dnew <- as.matrix(dist(xnew))
      anew <- sum(wmat * dnew)
      snew <- sum(wmat * (delta - dnew) ^ 2) / sum(wmat * (dnew - anew) ^ 2)
      if (verbose) {
        cat(
          "itel ",
          formatC(itel, format = "d"),
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
      sold <- snew
      dold <- dnew
      xold <- xnew
      aold <- anew
      itel <- itel + 1
    }
    return(list(
      x = xnew,
      s = snew,
      d = dnew,
      b = bmat,
      m = mmat,
      w = wmat,
      a = anew,
      u = umat,
      itel = itel
    ))
  }

torgerson <- function(delta, ndim) {
  dd <- delta ^ 2
  rd <- apply(dd, 1, mean)
  rr <- mean(dd)
  cc <- -.5 * (dd - outer(rd, rd, "+") + rr)
  ec <- eigen(cc)
  xx <- ec$vectors[, 1:ndim] %*% diag(sqrt(ec$values[1:ndim]))
  return(xx)
}
