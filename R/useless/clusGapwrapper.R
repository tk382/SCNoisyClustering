clusGapwrapper = function(X, numClust){
  tmpX        = X
  nClust      = numClust
  allk        = seq(3, 60, by=3); sigma = seq(0.05, 0.3, by=0.05);
  P           = multiple_kernel_tae(t(tmpX), allk, sigma)
  n           = nuclear_objective_c(P, 1)
  t           = tsne_spectral(n$S, numClust)
  return(list(cluster = t))
}
