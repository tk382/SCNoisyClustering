makeK = function(X, allk, sigma, bestind){
  #make multi-view kernels
  skip = round(nrow(X)/3)
  X1 = X[1:skip, ]
  X2 = X[(skip+1) : (2*skip),]
  X3 = X[(2*skip+1) : nrow(X), ]
  K1 = multiple_kernel_new(t(X1), allk, sigma)[[bestind]]
  K2 = multiple_kernel_new(t(X2), allk, sigma)[[bestind]]
  K3 = multiple_kernel_new(t(X3), allk, sigma)[[bestind]]
  K = array(0, dim=c(ncol(X), ncol(X), 3))
  K[,,1] = as.matrix(K1)
  K[,,2] = as.matrix(K2)
  K[,,3] = as.matrix(K3)
  rm(K1, K2, K3)
  #T is similarity matrices created from different features
  T = array(0, dim = c(ncol(X), ncol(X), 3))
  for (i in 1:3){
    Ki = K[,,i]
    D = diag(1/rowSums(Ki))
    Ti = D %*% Ki
    T[,,i] = Ti
  }
  D = diag(1/rowSums(Ki))
  Ti = D %*% Ki
  T[,,i] = Ti
  return(T)
}
