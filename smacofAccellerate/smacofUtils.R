
mPrint <- function(x,
                   digits = 10,
                   width = 15,
                   format = "f",
                   flag = "+") {
  print(noquote(
    formatC(
      x,
      digits = digits,
      width = width,
      format = format,
      flag = flag
    )
  ))
}

smacofLinePrint <- function(itel, sold, snew, cnew, labd, wd, dg) {
  cat(
    "itel",
    formatC(itel, width = wd, format = "d"),
    "sold",
    formatC(sold, digits = dg, format = "f"),
    "snew",
    formatC(snew, digits = dg, format = "f"),
    "chng",
    formatC(cnew, digits =  dg, format = "f"),
    "labd",
    formatC(labd, digits =  dg, format = "f"),
    "\n"
  )
}

smacofMakeBasis <- function(n, ndim, wgth = 1 - diag(n)) {
  vmat <- -wgth
  diag(vmat) <- -rowSums(vmat)
  y <- as.list(1:ndim)
  for (s in 0:(ndim - 1)) {
    ns <- n - s
    aux <- qr.Q(qr(ns * diag(ns) - 1))[, -ns]
    aux <- rbind(matrix(0, s, ns - 1), aux)
    sux <- crossprod(aux, vmat %*% aux)
    y[[s + 1]] <- tcrossprod(aux, solve(chol(sux)))
  }
  return(y)
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