data(ekman, package = "smacof")
delta <- as.matrix(1 - ekman)
wgth <- 1 / ((delta ^ 2) + diag(14))
diag(wgth) <- 0
vmat <- -wgth
diag(vmat) <- -rowSums(vmat)
h1<-smacofAccelerate(delta, wgth =  wgth, verbose = 0, itmax = 10000)
h1jacob <- smacofJacobian(h1$x, delta, wgth)
h1nacob <- numHess(h1$x, delta, wgth = wgth)
h1evecs <- eigen(h1jacob)$vectors
xx <- as.vector(vmat %*% h1$x)
print(crossprod(h1evecs, xx))
