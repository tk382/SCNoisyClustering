normalized_ng = function(S, numClust){
  rs = rowSums(S)
  L = diag(1/sqrt(rs)) %*% (diag(rs)- S) %*% diag(1/sqrt(rs))
  tmp = svd(L)
  U = tmp$v[, ((ncol(L)-numClust+1) : ncol(L))]
  norm_mat = matrix(rep(sqrt(rowSums(U^2)),numClust), ncol=numClust)
  errorind = which(norm_mat[,1]==0)
  norm_mat[errorind, ] = 1
  U = U/norm_mat
  C = kmeans(U, numClust, nstart = 200, iter.max=20)
  return(C$cluster)
}
