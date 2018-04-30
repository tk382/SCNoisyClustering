#tsne_spectral = function(A, L, numClust){
tsne_spectral = function(A, numClust, numEigen = NA){
  # D = diag(apply(S,MARGIN=2,FUN=sum))
  # L = D - S
  # compute the eigenvalues and eigenvectors of P
  if(is.na(numEigen)){numEigen = numClust}
  rs = rowSums(A)
  L = diag(1/sqrt(rs)) %*% (diag(rs) - A) %*% diag(1/sqrt(rs))
  eigen_L = eigen(L)
  U = eigen_L$vectors
  D = eigen_L$values
  U_index = seq(ncol(U),(ncol(U) - numEigen + 1))
  F_last = tsne(A, k = numEigen, initial_config = U[,U_index])
  C = kmeans(F_last, numClust, nstart=200)
  return(C$cluster)
}
