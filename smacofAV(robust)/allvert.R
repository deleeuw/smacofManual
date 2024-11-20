perms <- matrix(c(1,2,3,1,3,2,2,1,3,2,3,1,3,1,2,3,2,1), 6, 3, byrow = TRUE)
delta <- c(1, 4, 2)
dbase <- function(x, delta) {
  a <- matrix(0, 3, 2)
  a[x[1], ] <- c(0, 0)
  a[x[2], ] <- c(0, 1)
  a[x[3], ] <- c(1, 1)
  d <- matrix(0, 3, 2)
  d[1, ] <- abs(a[1, ] - a[2,])
  d[2, ] <- abs(a[1, ] - a[3,])
  d[3, ] <- abs(a[2, ] - a[3,])
  d <- rbind(diag(2), d)
  e <- c(0, 0, delta)
  for (i in 1:4) {
    for
  }
 }
