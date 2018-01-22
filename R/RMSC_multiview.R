RMSC_multiview = function(X, lambda, numClust, mu=1e-3, rho = 1.9, max_iter=100, eps=1e-9, verbose=FALSE){
  if(verbose){
    print('splitting data and computing kernel..')
  }
  skip = round(nrow(X)/3)
  X1 = X[1:skip, ]
  X2 = X[(skip+1) : (2*skip),]
  X3 = X[(2*skip+1) : nrow(X), ]
  K1 = multiple_kernel_new(t(X1), 1)[[21]]
  K2 = multiple_kernel_new(t(X2), 1)[[21]]
  K3 = multiple_kernel_new(t(X3), 1)[[21]]
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

  dims = dim(T);
  m = dims[1]; p = dims[2]; n = dims[3];
  if(m!=p){
    stop('input matrix T must be square transition matrix')
  }
  Z     = array(0,dim=c(m,p))
  E     = array(rnorm(m*p*n), dim = c(m,p,n))
  Y     = array(0, dim=c(m,p,n))
  B     = array(0, dim=c(m,p,n))
  Q     = array(0, dim=c(m,p))
  P     = array(0, dim=c(m,p))
  P_old = matrix(rnorm(m*p), m, p)
  L     = eigen((P_old+t(P_old))/2)$vectors[, 1:(numClust)]
  e     = rep(1, m)

  step = 0
  while(1){
    step = step + 1
    max_inf_norm = -1
    for (i in 1:n){
      Ti = T[,,i]
      Ei = E[,,i]
      diff = Ti-Ei-P
      inf_norm = norm(diff, 'I')
      max_inf_norm = max(max_inf_norm, inf_norm)
    }
    funV = sum(svd(P)$d) + lambda*sum(abs(as.numeric(E)))
    relChg = norm(P-P_old, 'F')/max(1, norm(P_old, 'F'))
    P_old = P;
    tmp = P-Q
    max_inf_norm2 = norm(tmp, 'I')
    if(verbose){
      cat(paste0('iter',step,': \n max_inf_norm=',round(max_inf_norm,4),
            ',\n relChg=',round(relChg,4), ',\n mu=',round(mu,4),', \n inf_norm2=', round(max_inf_norm2,4),
            ',\n funV=',round(funV,4), '\n'))
    }
    if (step > 1 && max_inf_norm < eps){
      break;
    }
    if (step > max_iter){
      print(paste('reached max iteration : ', step))
      break;
    }

    #update P
    B = 1/(n+1) * (Q-Z/mu + apply(T-E-Y/mu, c(1,2),sum))
    P = nonnegASC(B)
    pi = irlba(t(P), 1)$v
    Dist = as.numeric(pi)/sum(pi)
    Dist2 = sum(pi)/as.numeric(pi)
    sspi = sqrt(diag(Dist2))
    spi = sqrt(diag(Dist))
    P = (spi %*% P %*% sspi + sspi %*% t(P) %*% spi)/2

    #update Q
    Q = updateQ(P,Z,mu,m)

    #update Ei and Yi
    for (i in 1:n){
      C = T[,,i]-P-Y[,,i]/mu
      E[,,i] = pmax(C-lambda/mu,0) + pmin(C+lambda/mu,0)
      Y[,,i] = Y[,,i] + mu * (P + E[,,i] - T[,,i]);
    }

    #update Z
    Z = Z+mu*(P-Q)

    #update mu
    mu = min(rho*mu, 1e+10)
  }

  pi = irlba(t(P), 1)$v
  Dist = as.numeric(pi)/sum(pi)
  Dist2 = sum(pi)/as.numeric(pi)
  sspi = sqrt(diag(Dist2))
  spi = sqrt(diag(Dist))
  P = (spi %*% P %*% sspi + sspi %*% t(P) %*% spi)/2
  return(list(P=P, E=E))
}

updateQ = function(P,Z,mu,m){
  M = P+Z/mu
  C = 1/mu
  USigV = svd(M)
  U = USigV$u; Sigma = USigV$d; V = USigV$v;
  svp = sum(Sigma > C)
  if(svp>=2){
    Sigma = Sigma[1:svp]-C;
    Q = U[, 1:svp] %*% diag(Sigma) %*% t(V[,1:svp])
  }else if(svp==1){
    Sigma = Sigma[1]-C;
    Q = as.matrix(U[,1]) %*% matrix(Sigma, nrow=1) %*%  as.matrix(t(V[,1]))
  }else{
    svp=1
    Q = matrix(0,m,m)
  }
  return(Q)
}

# newupdateQ = function(P,Z,L,mu,m){
#   M = t(L) %*% (P + Z/mu) %*% L
#   C = 1/mu
#   USigV = svd(M)
#   U = USigV$u; Sigma = USigV$d; V = USigV$v;
#   svp = sum(Sigma > C)
#   if(svp>=2){
#     Sigma = Sigma[1:svp]-C;
#     Q = U[, 1:svp] %*% diag(Sigma) %*% t(V[,1:svp])
#   }else if(svp==1){
#     Sigma = Sigma[1]-C;
#     Q = as.matrix(U[,1]) %*% matrix(Sigma, nrow=1) %*%  as.matrix(t(V[,1]))
#   }else{
#     svp=1
#     Q = matrix(0,ncol(L),ncol(L))
#   }
#   Q = L %*% Q %*% t(L)
#   return(Q)
# }
