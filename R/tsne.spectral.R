#' dimension reduction wrapper
#'
#' @param A Estimated similarity matrix.
#' @param numClust Number of clusters. This can't be NA. If unknown, it should be estimated from getNumClust function using A
#' @param numEigen number of eigenvectors to be used. If NA, the number of clusters will be used
#' @examples
#'
#' @export
tsne.spectral = function(A, numClust, numEigen = NA){
  # compute the eigenvalues and eigenvectors of P
  if(is.na(numEigen)){numEigen = numClust}
  rs = rowSums(A) + 1
  L = diag(1/sqrt(rs)) %*% (diag(rs) - A) %*% diag(1/sqrt(rs))
  L = (L + t(L))/2
  eigen_L = eigen(L)
  U = eigen_L$vectors
  U_index = seq(ncol(U), (ncol(U)-numEigen+1))
  F_last = tsne_c(A, initial_config = U[,U_index], k = numEigen)
  C = kmeans(F_last, numClust, nstart=200)
  return(list(cluster = C$cluster, tsne = F_last))
}

