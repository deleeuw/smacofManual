huber <- function(x, c) {
  if (abs(x) < c) {
    return((x ^ 2) / 2)
  }
  else {
    return(c * abs(x) - (c ^ 2) / 2)
  }
}

huber_deriv <- function(x, c) {
  if (abs(x) < c) {
    return(x)
  }
  else {
    return(sign(x) * c)
  }
}

huber_majorize <- function(x, y, c) {
  return((huber_deriv(y, c) / (2 * y)) * (x ^ 2 - y ^ 2) + huber(y, c))
}

plotHuber <- function(y, c) {
  x <- seq(-3, 3, length = 100)
  z <- rep(0, 100)
  for (i in 1:100) z[i] <- huber(x[i], c)
  plot(x, z, type = "l", col = "RED", lwd = 2)
  for (i in 1:100) z[i] <- huber_majorize(x[i], y, c)
  lines(x, z, col = "BLUE", lwd = 2)
  abline(v = c)
  abline(v = -c)
  abline(v = y, col = "GREEN", lwd = 3)
}