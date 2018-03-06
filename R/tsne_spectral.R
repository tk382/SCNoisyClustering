tsne_spectral = function(S, numClust){
  D = diag(apply(S,MARGIN=2,FUN=sum))
  L = D - S
  # compute the eigenvalues and eigenvectors of P
  eigen_L = eigen(L)
  U = eigen_L$vectors
  D = eigen_L$values
  U_index = seq(ncol(U),(ncol(U) - numClust + 1))
  F_last = tsne(S, k = numClust, initial_config = U[,U_index])
  C = kmeans(F_last, numClust, nstart=200)
  return(C$cluster)
}
