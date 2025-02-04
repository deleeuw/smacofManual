library(RSpectra)

rderiv1 <- function(r, de = delta, w = wght, dd = dx) {
  n <- nrow(delta)
  dr <- de ^ r
  lr <- log(de + diag(n))
  return(sum(w * dr * lr * (dr - dd))) 
}

rderiv2 <- function(r, de = delta, w = wght, dd = dx) {
  n <- nrow(delta)
  dr <- de ^ r
  lr <- log(de + diag(n)) ^ 2
  return(sum(w * dr * lr * (2 * dr - dd)))
}

rfunc <- function(r, de = delta, w = wght, dd = dx) {
  n <- nrow(de)
  dr <- de ^ r
  return(sum(w * (dr - dd) ^ 2))
}

torgerson <- function(delta, ndim = 2) {
  dd <- delta ^ 2
  rd <- apply(dd, 1, mean)
  md <- mean(dd)
  ed <- -.5 * (dd - outer(rd, rd, "+") + md)
  ee <- eigs_sym(ed, ndim)
  return(ee$vectors %*% diag(sqrt(ee$values)))
}

data(ekman, package = "smacof")
delta <- as.matrix(1 - ekman)
dx <- as.matrix(dist(torgerson(delta)))
wght <- 1 - diag(14)