n <- 4
d <- 1:6
set.seed(12345)
x <- matrix(rnorm(8), 4, 2)
s <- tcrossprod(x)
v <- as.vector(s)

smacofMakeAij <- function(i, j, n) {
  ei <- ifelse(i == 1:n, 1, 0)
  ej <- ifelse(j == 1:n, 1, 0)
  return(outer(ei - ej, ei - ej))
}

a <- matrix(0, n * (n - 1) / 2, n ^ 2)
k <- 1
for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
    a[k, ] <- as.vector(smacofMakeAij(i, j, n))
    k <- k + 1
  }
}

kk <- matrix(0, n ^ 2, n ^ 2)
for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
    aij <- smacofMakeAij(i, j, n)
    kk <- kk + kronecker(aij, aij)
  }
}

first <- sum(v * (kk %*% v))

second <- third <- 0
for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
    aij <- smacofMakeAij(i, j, n)
    second <- second + sum(diag(aij %*% s %*% aij %*% s))
    third <- third + sum(aij * s) ^ 2
  }
}

cs <- qr.solve(a, d)

cc  <- matrix(cs, 4, 4)

cc <- (cc + t(cc)) / 2
