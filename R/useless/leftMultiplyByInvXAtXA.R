leftMultiplyByInvXAtXA <- function (n, ind, val, phi = rep(1, n), w = defaultWeights(n), verbose = FALSE)
{
  a <- dim(val)[1]
  p <- dim(val)[2]
  o <- order(ind)
  ind <- ind[o]
  val <- val[o, , drop = FALSE]
  r <- matrix(numeric(a * p), nrow = a, ncol = p)
  if (length(w) != (n - 1)) {
    stop("'w' needs to be of length n-1")
  }
  if (a != 0) {
    v <- diff(c(0, cumsum(phi^2)[ind], n))
    d <- w[ind]
    R <- matrix(numeric((a + 2) * p), ncol = p)
    val <- apply(val, 2, function(x) {
      x/d
    })
    R[1, ] <- numeric(p)
    R[2:(a + 1), ] <- val
    R[(a + 2), ] <- numeric(p)
    gamma <- apply(R, 2, diff)
    delta <- gamma/v
    r <- -diff(delta)
    r <- r/d
  }
  return(r)
}
