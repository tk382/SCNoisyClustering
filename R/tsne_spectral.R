#tsne_spectral = function(A, L, numClust){
tsne_spectral = function(A, numClust, numEigen = NA){
  # compute the eigenvalues and eigenvectors of P
  if(is.na(numEigen)){numEigen = numClust}
  rs = rowSums(A) + 1
  L = diag(1/sqrt(rs)) %*% (diag(rs) - A) %*% diag(1/sqrt(rs))
  L = (L + t(L))/2
  eigen_L = eigen(L)
  U = eigen_L$vectors
  #D = eigen_L$values
  U_index = seq(ncol(U), (ncol(U)-numEigen+1))
  F_last = tsne_c(A, initial_config = U[,U_index], k = numEigen)
  C = kmeans(F_last, numClust, nstart=200)
  return(list(cluster = C$cluster, tsne = F_last))
}

