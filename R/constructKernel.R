constructKernel = function(X, sigma){
  out = exp(-as.matrix(dist(X)) / (2 * sigma^2))
  return(out)
}
