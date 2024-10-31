phfunc <- function() {
  plot(1, type = "n", xlim = c(-5, 5), ylim = c(0, 3), xlab = "x", ylab = "h")
  lines(matrix(c(-5, 0, -1, 0), 2, 2, byrow = TRUE), col = "RED", lwd = 2)
  lines(matrix(c(-1, 1, 0, 1), 2, 2, byrow = TRUE), col = "RED", lwd = 2)
  lines(matrix(c(0, 3, 2, 3), 2, 2, byrow = TRUE), col = "RED", lwd = 2)
  lines(matrix(c(2, 2, 4, 2), 2, 2, byrow = TRUE), col = "RED", lwd = 2)
  lines(matrix(c(4, 0, 5, 0), 2, 2, byrow = TRUE), col = "RED", lwd = 2)
  return()
}

pgfunc <- function() {
  plot(1, type = "n", xlim = c(-5, 5), ylim = c(-2, 8), xlab = "x", ylab = "h")
  lines(matrix(c(-5, 0, -1, 0), 2, 2, byrow = TRUE), col = "RED", lwd = 2)
  lines(matrix(c(-1, -1, 0, 0), 2, 2, byrow = TRUE), col = "RED", lwd = 2)
  lines(matrix(c(0, 0, 2, 6), 2, 2, byrow = TRUE), col = "RED", lwd = 2)
  lines(matrix(c(2, 4, 4, 8), 2, 2, byrow = TRUE), col = "RED", lwd = 2)
  lines(matrix(c(4, 0, 5, 0), 2, 2, byrow = TRUE), col = "RED", lwd = 2)
  return()
}

pffunc <- function() {
  plot(1, type = "n", xlim = c(-5, 5), ylim = c(-2, 20), xlab = "x", ylab = "h")
  lines(matrix(c(-6, 0, -1, 0), 2, 2, byrow = TRUE), col = "RED", lwd = 2)
  x <- seq(-1, 0, length = 100)
  y <- (x ^ 2 - 1) / 2
  lines(x, y, col = "RED", lwd = 2)
  x <- seq(0, 2, length = 100)
  y <- (3 * x ^ 2 - 1) / 2
  lines(x, y, col = "RED", lwd = 2)
  x <- seq(2, 4, length = 100)
  y <- x ^ 2 + 1.5
  lines(x, y, col = "RED", lwd = 2)
  lines(matrix(c(4, 17.5, 6, 17.5), 2, 2, byrow = TRUE), col = "RED", lwd = 2)
  return()
}

x <- seq(-2 * pi, 2 * pi, length = 1000)
f <- g <- h <- rep(0, 1000)
for (i in 1:1000) {
  z <- x[i]
  f[i] <- ffunc(z)
  g[i] <- gfunc(z)
  h[i] <- hfunc(z)
}

fmaj <- function(y) {
  a <- ffunc(y) + hfunc(y) * (x ^ 2 - y ^ 2) / 2
  return(a)
}

ffunc <- function(z) {
  if (z < -1) {
    return(0.0)
  }
  if (z < 0.0) {
    return((z ^ 2 - 1) / 2.0)
  }
  if (z < 2.0) {
    return((3 * z ^ 2 - 1) / 2.0)
  }
  if (z < 4.0) {
    return(z ^ 2 + 1.5)
  }
  return(17.5)
}

gfunc <- function(z) {
  if (z < -1) {
    return(0.0)
  }
  if (z < 0.0) {
    return(z)
  }
  if (z < 2.0) {
    return(3.0 * z)
  }
  if (z < 4.0) {
    return(2.0 * z)
  }
  return(0.0)
}

hfunc <- function(z) {
  if (z < -1) {
    return(0.0)
  }
  if (z < 0.0) {
    return(1.0)
  }
  if (z < 2.0) {
    return(3.0)
  }
  if (z < 4.0) {
    return(2.0)
  }
  return(0.0)
}