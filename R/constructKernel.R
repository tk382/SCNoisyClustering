constructKernel = function(X, sigma){
  n = ncol(X)
  out = exp(-as.matrix(dist(X)) / (2 * sigma^2))
  return(out)
}
