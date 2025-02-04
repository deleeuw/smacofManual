library(RSpectra)

smacofPO <-
  function(delta,
           wght = 1 - diag(nrow(delta)),
           interval = c(.5, 2),
           fixed = FALSE,
           xe = NULL,
           ndim = 2,
           method = "newton",
           itmax = 1000,
           eps = 1e-10,
           verbose = TRUE) {
    nobj <- nrow(delta)
    vmat <- -wght
    diag(vmat) <- -rowSums(vmat)
    vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
    if (length(interval) == 1) {
      rold <- interval
    } else {
      rold <- (interval[1] + interval[2]) / 2
    }
    if (is.null(xe)) {
      xe <- smacofTorgerson(delta ^ rold, ndim)
    }
    dx <- as.matrix(dist(xe))
    g <- function(r, delta, de) {
      return(sum(((delta^r) - de)^2))
    }
    sold <- rfunc(rold, delta, wght, dx)
    itel <- 1
    repeat {
      e <- (delta + diag(nobj))^rold
      diag(e) <- 0
      b <- -wght * e * ifelse(dx == 0, 0, 1 / dx)
      diag(b) <- -rowSums(b)
      xe <- vinv %*% (b %*% xe)
      dx <- as.matrix(dist(xe))
      smid <- rfunc(rold, delta, wght, dx)
      if ((!fixed) && (method == "newton")) {
        rnew <- rNewton(
          rold,
          delta,
          wght,
          dx,
          itmax = 10,
          eps = 1e-10,
          verbose = FALSE
        )$r
      }
      if ((!fixed) && (method == "optimize")) {
        rnew <- optimize(g,
                         interval = interval,
                         delta = delta,
                         de = dx)$minimum
      }
      if (fixed) {
        rnew <- rold
      }
      snew <- rfunc(rnew, delta, wght, dx)
      if (verbose) {
        cat(
          "itel ",
          formatC(itel, format = "d"),
          "sold ",
          formatC(sold, digits = 10, format = "f"),
          "smid ",
          formatC(smid, digits = 10, format = "f"),
          "snew ",
          formatC(snew, digits = 10, format = "f"),
          "\n"
        )
      }
      if (((sold - snew) < 1e-10) || (itel == itmax)) {
        break
      }
      itel <- itel + 1
      rold <- rnew
      sold <- snew
    }
    return(list(
      x = xe,
      d = dx,
      e = delta^rnew,
      r = rnew,
      itel = itel,
      stress = snew
    ))
  }

rderiv1 <- function(r,
                    de = delta,
                    w = wght,
                    dd = dx) {
  n <- nrow(de)
  dr <- (de + diag(n))^r
  diag(dr) <- 0
  lr <- log(de + diag(n))
  return(sum(w * dr * lr * (dr - dd)))
}

rderiv2 <- function(r,
                    de = delta,
                    w = wght,
                    dd = dx) {
  n <- nrow(de)
  dr <- (de + diag(n))^r
  diag(dr) <- 0
  lr <- log(de + diag(n))^2
  return(sum(w * dr * lr * (2 * dr - dd)))
}

rfunc <- function(r,
                  de = delta,
                  w = wght,
                  dd = dx) {
  n <- nrow(de)
  dr <- (de + diag(n))^r
  diag(dr) <- 0
  return(sum(w * (dr - dd)^2))
}

rNewton <- function(rold,
                    delta,
                    wght,
                    dx,
                    itmax = 10,
                    eps = 1e-10,
                    verbose = TRUE) {
  itel <- 1
  fold <- rfunc(rold, de = delta, w = wght, dd = dx)
  repeat {
    g <- rderiv1(rold,
                 de = delta,
                 w = wght,
                 dd = dx)
    h <- rderiv2(rold,
                 de = delta,
                 w = wght,
                 dd = dx)
    u <- g / h
    rnew <- rold - u
    fnew <- rfunc(rnew,
                  de = delta,
                  w = wght,
                  dd = dx)
    if (verbose) {
      print(c(rold, rnew, fold, fnew), digits = 10)
    }
    if ((itel == itmax) || (abs(g) < eps)) {
      break
    }
    itel <- itel + 1
    fold <- fnew
    rold <- rnew
  }
  return(list(r = rnew, f = fnew, itel = itel))
}

smacofTorgerson <- function(delta, ndim = 2) {
  dd <- delta^2
  rd <- apply(dd, 1, mean)
  md <- mean(dd)
  ed <- -.5 * (dd - outer(rd, rd, "+") + md)
  ee <- eigs_sym(ed, ndim)
  return(ee$vectors %*% diag(sqrt(ee$values)))
}

smacofLTR <- function(x) {
  n <- nrow(x)
  return(x[outer(1:n, 1:n, "<")])
}

smacofPOShepardPlot <- function(delta, d, r) {
  dt <- smacofLTR(delta)
  dd <- smacofLTR(d)
  plot(dt,
       dd,
       xlab = "delta",
       ylab = "dist",
       col = "BLUE")
  ds <- sort(dt)
  lines(ds, ds^r, col = "RED", lwd = 3)
}


smacofPORplot <- function(delta,
                          wght = 1 - diag(nrow(delta)),
                          rmin = 1,
                          rmax = 2,
                          rlength = 100,
                          ndim = 2,
                          itmax = 1000,
                          eps = 1e-10,
                          verbose = FALSE
                          ) {
  r <- seq(rmin, rmax, length = rlength)
  s <- rep(0, rlength)
  for (i in 1:rlength) {
    if (i == 1) {
      xe <- NULL
    }
    h <- smacofPO(
      delta,
      xe = xe,
      interval = r[i],
      fixed = TRUE,
      eps = eps,
      ndim = ndim,
      itmax = itmax,
      verbose = verbose
    )
    s[i] <- h$stress
    xe <- h$x
    if (i == 1) {
      plot(xe, col = "RED")
    } else {
      points(xe)
    }
  }
  plot(
    r,
    s,
    type = "l",
    lwd = 2,
    col = "RED",
    ylab = "stress"
  )
}
