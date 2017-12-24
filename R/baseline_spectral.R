baseline_spectral = function(Phat, numClust){
  tmp = irlba(Phat, numClust)
  U = tmp$u
  norm_mat = matrix(rep(sqrt(rowSums(U^2)),numClust), ncol=numClust)
  errorind = which(norm_mat[,1]==0)
  norm_mat[errorind, ] = 1
  U = U/norm_mat
  C = kmeans(U, numClust)
  return(C$cluster)
}
