## This function minimizes f(x) = c - b'x + x'Ax
## over x ≥ 0. A is assumed to be pds.
# It needs an initial estimate for x and an upper bound for the
# largest eigenvalue of A.



smacofNonnegativeLeastSquares <- function(a,
                                          b,
                                          c,
                                          xold = rep(0, nrow(a)),
                                          bnd = max(rowSums(abs(a))),
                                          itmax = 100,
                                          eps = 1e-10,
                                          verbose = TRUE) {
  itel <- 1
  fold <- c - sum(b * xold) + sum(xold * (a %*% xold))
  repeat {
    z <- xold - (a %*% xold - b) / bnd
    xnew <- pmax(0, z)
    fnew <- c - sum(b * xnew) + sum(xold * (a %*% xnew))
    if (verbose) {
      cat(
        "itel ",
        formatC(itel, format = "d", width = 4),
        "fold ",
        formatC(fold, format = "f", digits = 10),
        "fnew ",
        formatC(fnew, format = "f", digits = 10),
        "\n"
      )
    }
    if ((itel == itmax) || ((fold - fnew) < eps)) {
      break
    }
    itel <- itel + 1
    xold <- xnew
    fold <- fnew
  }
}