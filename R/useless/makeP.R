makeP = function(X, allk, sigma){
  K = multiple_kernel_new(t(X), allk, sigma, cores.ratio=1)
  P = array(0,dim=c(ncol(X),ncol(X),length(K)))
  for (i in 1:length(K)){
    P[,,i] = as.matrix(K[[i]])
  }
  return(P)
}
