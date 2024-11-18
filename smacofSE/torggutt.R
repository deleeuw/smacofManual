gruijter <-
  structure(
    c(
      5.63,
      5.27,
      4.6,
      4.8,
      7.54,
      6.73,
      7.18,
      6.17,
      6.72,
      5.64,
      6.22,
      5.12,
      4.59,
      7.22,
      5.47,
      5.46,
      4.97,
      8.13,
      7.55,
      6.9,
      4.67,
      3.2,
      7.84,
      6.73,
      7.28,
      6.13,
      7.8,
      7.08,
      6.96,
      6.04,
      4.08,
      6.34,
      7.42,
      6.88,
      6.36,
      7.36
    ),
    Labels = c("KVP", "PvdA", "VVD",
               "ARP", "CHU", "CPN", "PSP", "BP", "D66"),
    Size = 9L,
    call = quote(as.dist.default(m = polpar)),
    class = "dist",
    Diag = FALSE,
    Upper = FALSE
  )

labels <- c("KVP", "PvdA", "VVD",
            "ARP", "CHU", "CPN", "PSP", "BP", "D66")

gruijter <- as.matrix(gruijter)

gruijter[3, 2] <- gruijter[2, 3] <- NA
gruijter[3, 5] <- gruijter[5, 3] <- NA
gruijter[7, 4] <- gruijter[4, 7] <- NA
gruijter[2, 9] <- gruijter[9, 2] <- NA


n <- 9
tlist <- list()
nmis <- 0
for (j in 1:(n - 1)) {
  for (i in (j + 1):n) {
    if (is.na(gruijter[i, j])) {
      e <- matrix(0, n, n)
      e[i, j] <- e[j, i] <- 1
      tlist <- c(tlist, list(-smacofDoubleCenter(e) / 2))
      nmis <- nmis + 1
    }
  }
}
cc <- matrix(0, nmis, nmis)
for (i in 1:nmis) {
  for (j in 1:nmis) {
    cc[i, j] <- sum(tlist[[i]] * tlist[[j]])
  }
}
  
