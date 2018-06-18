leftMultiplyByXt <- function (Y, phi = rep(1, nrow(Y)), w = defaultWeights(nrow(Y)), verbose = FALSE)
{
  n <- as.numeric(dim(Y)[1])
  p <- ncol(Y)

  Y <- Y * phi

  u <- apply(Y, 2, cumsum)
  if (length(w) != (n - 1)) {
    stop("w needs to be of length nrow(Y)-1")
  }
  C <- apply(u, 2, function(x) {
    w * (cumsum(phi^2)[1:(n - 1)] * x[n]/n - x[1:(n - 1)])
  })
  dim(C) <- c(n - 1, p)
  return(C)
}
