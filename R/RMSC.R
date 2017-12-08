library(matrixStats)
library(irlba)
T1 = as.matrix(read.csv('../RMSC/T1.csv', header = FALSE))
T2 = as.matrix(read.csv('../RMSC/T2.csv', header = FALSE))
T3 = as.matrix(read.csv('../RMSC/T3.csv', header = FALSE))
T = array(0,dim=c(2000,2000,3))
T[,,1] = T1; T[,,2] = T2; T[,,3] = T3
lambda = 0.005
verbose=TRUE
eps = 1e-6
max_iter = 300
rm(T1, T2, T3)

RMSC = function(T, lambda, mu=1e-3, rho = 1.9, max_iter=100, eps=1e-9, verbose=FALSE){
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
  e     = rep(1, m)

  step = 0
  funVal = numeric()
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
      print(paste('iter',step,': max_inf_norm=',max_inf_norm,
            ', relChg=',relChg, ',mu=',mu,', inf_norm2=', max_inf_norm2,
            ', funV=',funV))
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
    for (i in 1:m){
      if (sum(P[i,])-1.0 > 1e-10){stop('sum to 1 error')}
    }

    #update Q
    M = P+Z/mu
    C = 1/mu
    USigV = svd(M)
    U = USigV$u; Sigma = USigV$d; V = USigV$v;
    svp = sum(Sigma > C)
    if(svp>=2){
      Sigma = Sigma[1:svp]-C;
      Q = U[, 1:svp] %*% diag(Sigma) %*% t(V[,1:svp])
    }else if(svp==1){
      svp=1
      Sigma = Sigma[1]-C;
      Q = as.matrix(U[,1]) %*% matrix(Sigma, nrow=1) %*%  as.matrix(t(V[,1]))
    }else{
      svp=1
      Sigma = 0
      Q = matrix(0,m,m)
    }
    #update Ei
    for (i in 1:n){
      C = T[,,i]-P-Y[,,i]/mu
      E[,,i] = pmax(C-lambda/mu,0)+pmin(C+lambda/mu,0)
      Y[,,i] = Y[,,i] + mu*(P+E[,,i]-T[,,i]);
    }
    Z = Z+mu*(P-Q)
    mu = min(rho*mu, 1e+10)
  }
  pi = irlba(t(P), 1)$v
  Dist = pi / sum(pi)
  pi = (diag(as.numeric(Dist)))
  sspi = sqrt(solve(pi))
  spi = sqrt(pi)
  P = (spi %*% P %*% sspi + sspi %*% t(P) %*% spi)/2
  return(P)
}
