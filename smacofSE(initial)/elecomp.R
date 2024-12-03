smacofMakeAij <- function(i, j, n) {
  ei <- ifelse(i == 1:n, 1, 0)
  ej <- ifelse(j == 1:n, 1, 0)
  return(outer(ei - ej, ei - ej))
}
n1 <- 20
n2 <- n1 ^ 2
n3 <- n1 * (n1 - 1) / 2
w <- 1 - diag(n1)
set.seed(12345)
w <- matrix(sample(1:2, n2, replace = TRUE), n1, n1)
w <- floor((w + t(w)) / 2)
diag(w) <- 0
aa <- matrix(0, n2, n3)
v1 <- matrix(0, n2, n2)
v2 <- matrix(0, n2, n2)
w2 <- rep(0, n3)
k <- 1
for (j in 1:n1) {
  for (i in 1:n1) {
    if (i >= j) {
      next
    }
    aij <- smacofMakeAij(i, j, n1)
    aa[, k] <- as.vector(aij)
    w2[k] <- sqrt(w[i, j])
    v1 <- v1 + w[i, j] * outer(aa[, k],aa[, k])
    v2 <- v2 + w[i, j] * kronecker(aij, aij)
    k <- k + 1
  }
}
v3 <- matrix(0, n2, n2)
for (i in 1:n1) {
  indi <- 1:n1  + n1 * (i - 1)
  for (j in 1:n1) {
    if (i == j) {
      next
    }
    jndi <- 1:n1  + n1 * (j - 1)
    vs <- w[i, j] * smacofMakeAij(i, j, n1)
    v3[indi, jndi] <- -vs
    v3[indi, indi] <- v3[indi, indi] + vs
  }
}
ab <- aa %*% diag(w2)
v4 <- tcrossprod(ab)
uu <- crossprod(ab)

