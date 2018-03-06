baseline_spectral = function(S, numClust){
  #D = diag(apply(S,MARGIN=2,FUN=sum))
  #L = D - S
  tmp = irlba(S, numClust)
  U = tmp$u
  norm_mat = matrix(rep(sqrt(rowSums(U^2)),numClust), ncol=numClust)
  errorind = which(norm_mat[,1]==0)
  norm_mat[errorind, ] = 1
  U = U/norm_mat
  #C = kmeans(U, numClust, nstart = 200, iter.max=20)
  C = kmeans(U, numClust)
  return(C$cluster)
}

