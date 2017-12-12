make_T_matrix = function(X, shuffled, n){
  #shuffled is a index vector with random or user-specific (needs editing)
  skip = floor(nrow(X) / n)
  T = array(0, dim=c(ncol(X), ncol(X), n))
  for (i in 1:n){
    if(i!=n){
      Xi = X[shuffled[(skip*(i-1) + 1) : (skip * i)], ]
    }
    if(i==n){
      Xi = X[shuffled[(skip*(i-1) + 1) : nrow(X)], ]
    }
    Ki = constructKernel(t(Xi), 5)
    D = diag(1/rowSums(Ki))
    Ti = D %*% Ki
    T[,,i] = Ti
  }
  return(T)
}
