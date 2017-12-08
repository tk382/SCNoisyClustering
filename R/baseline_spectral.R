baseline_spectral = function(Phat, numClust, truth){
  tmp = irlba(Phat, numClust)
  U = tmp$u
  norm_mat = matrix(rep(sqrt(rowSums(U^2)),numClust), ncol=numClust)
  for (i in 1:nrow(norm_mat)){
    if(norm_mat[i,1]==0){norm_mat[i,] = 1}
  }
  U = U/norm_mat
  for (i in 1:10){
    C = kmeans(U, numClust)
  }
  return(C$cluster)
}
