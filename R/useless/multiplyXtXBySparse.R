multiplyXtXBySparse <- function (n, ind, val, phi = rep(1, n), wts = defaultWeights(n), verbose = FALSE)
{
  if (length(wts) != (n - 1)) {
    stop("Argument 'w' has to be of length nrow(Y)-1")
  }
  a <- nrow(val)
  p <- ncol(val)
  indrev <- (n-1):1
  if (a != 0) {
    o <- order(ind)
    ind <- ind[o]
    val <- val[o, , drop = FALSE]
    r <- val * wts[ind]
    s <- cumsum(phi^2)[ind] * r
    s <- colSums(s)/n
    T <- matrix(numeric((n - 1) * p), nrow = (n-1), ncol = p)
    T[indrev[ind], ] <- r[, , drop = FALSE]
    T <- apply(T, 2, cumsum)
    T <- T[indrev, , drop = FALSE]
    u <- sweep(T, 2, s)
    u <- u * (phi[-n])^2
    U <- apply(u, 2, cumsum)
    C <- apply(U, 2, function(x) x*wts)
  }
  return(C)
}
