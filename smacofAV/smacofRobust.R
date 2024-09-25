smacofRobust <- function(delta,
                         weights = 1 - diag(nrow(delta)),
                         engine = smacofAV,
                         ndim = 2,
                         cons = 0,
                         itmax = 1000,
                         eps = 1e-15,
                         verbose = TRUE) {
  nobj <- nrow(delta)
  wmax <- max(weights)
  xold <- smacofTorgerson(delta, ndim)
  dold <- as.matrix(dist(xold))
  h <- engine(nobj, weights, delta, dold, cons)
  rold <- h$resi
  sold <- sum(weights * rold)
  wold <- h$wght
  itel <- 1
  repeat {
    vmat <- -wold
    diag(vmat) <- -rowSums(vmat)
    vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
    bmat <- -wold * delta / (dold + diag(nobj))
    diag(bmat) <- -rowSums(bmat)
    xnew <- vinv %*% (bmat %*% xold)
    dnew <- as.matrix(dist(xnew))
    h <- engine(nobj, weights, delta, dnew, cons)
    rnew <- h$resi
    wnew <- h$wght
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

smacofTorgerson <- function(delta, ndim) {
  dd <- delta ^ 2
  rd <- apply(dd, 1, mean)
  md <- mean(dd)
  sd <- -.5 * (dd - outer(rd, rd, "+") + md)
  ed <- eigen(sd)
  return(ed$vectors[, 1:ndim] %*% diag(sqrt(ed$values[1:ndim])))
}

smacofAV <- function(nobj, wmat, delta, dmat, cons) {
  resi <- sqrt((delta - dmat) ^ 2 + cons)
  wght <- wmat / (resi + diag(nobj))
  return(list(resi = resi, wght = wght))
}

smacofConvolution <- function(nobj, wmat, delta, dmat, cons) {
  difi <- delta - dmat
  resi <- difi * (2 * pnorm(difi / cons) - 1) + 2 * cons * dnorm(difi / cons)
  wght <- wmat * (pnorm(difi / cons) - 0.5) / (difi + diag(nobj))
  return(list(resi = resi, wght = wght))
}

smacofTukey <- function(nobj, wmat, delta, dmat, cons) {
  difi <- delta - dmat
  resi <- ((cons ^ 2) / 6) * ifelse(abs(difi) < cons, (1 - (1 - (difi / cons) ^ 2) ^ 3), 1)
  wght <- ifelse(abs(difi) < cons, wmat * (1 - (difi / cons) ^ 2) ^ 2, 0) / 2
  return(list(resi = resi, wght = wght))
}

smacofHuber <- function(nobj, wmat, delta, dmat, cons) {
  difi <- delta - dmat
  resi <- ifelse(abs(difi) < cons, (difi ^ 2) / 2, cons * abs(difi) - ((cons ^ 2) / 2))
  wght <- ifelse(abs(difi) < cons, wmat,
                 wmat * sign(difi - cons) * cons / (difi + diag(nobj)))
  return(list(resi = resi, wght = wght))
}
