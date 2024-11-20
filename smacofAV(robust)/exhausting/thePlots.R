thePlots <- function(xmin = -3,
                     xmax = 3,
                     ymin = -3,
                     ymax = 3,
                     pow = 1,
                     n = 100) {
  x <- seq(xmin, xmax, length = n)
  y <- seq(ymin, ymax, length = n)
  z <- matrix(0, n, n)
  if (pow == 1) {
    for (i in 1:n) {
      for (j in 1:n) {
        z[i, j] <- z[i, j] + abs(1 - abs(y[j]))
        z[i, j] <- z[i, j] + abs(4 - abs(x[i] + y[j]))
        z[i, j] <- z[i, j] + abs(2 - abs(x[i]))
      }
    }
  }
  if (pow == 2) {
    for (i in 1:n) {
      for (j in 1:n) {
        z[i, j] <- z[i, j] + (1 - abs(y[j])) ^ 2
        z[i, j] <- z[i, j] + (4 - abs(x[i] + y[j])) ^ 2
        z[i, j] <- z[i, j] + (2 - abs(x[i])) ^ 2
      }
    }
  }
  return(list(x = x, y = y, z = z))
}
