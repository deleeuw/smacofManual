x <- c()
y <- c()
for (i in 1:1000) {
  h <- smacofRobust(delta, verbose = FALSE, cons = 0, ndim = 1, xold = matrix(rnorm(3), 3, 1))
  x <- c(x, h$s)
  y <- c(y, h$itel)
}

x <- y <- seq(-3, 3, length = 1000)
z <- matrix(0, 1000, 1000)
for (i in 1:1000) {
  for (j in 1:1000) {
    z[i,j] <- z[i, j] + abs(2 * y[j])
    z[i,j] <- z[i, j] + abs(y[j] - 3 * x[i])
    z[i,j] <- z[i, j] + abs(y[j] + 3 * x[i])
  }
}
