Laplacian1 = function(P, tau, max_iter=100, eps=1e-9, verbose=FALSE){
  # Observe unnormalized similarity matrices P[,,i]
  # First create Mi = Di^{-1/2}*(Di-Pi)*Di^{-1/2}
  # Then take the minimum of
  # tau*norm(I-L, nuclear) + \sum \|L-M\|_F^2 / (2n^2 \sigma_i^2) + \sigma_i / 2

  dims = dim(P)
  funV = rep(0,max_iter)
  m = dims[1]; p = dims[2]; n = dims[3];
  if(m!=p){
    stop('input matrix P must be square transition matrix')
  }
  #compute array M
  M = H = P
  for (i in 1:n){
    rs = pmax(rowSums(P[,,i]), 1)
    rsinv = 1/sqrt(rs)
    DD    = matrix(rep(rsinv, m), nrow=m); DD = DD*t(DD)
    M[,,i] = (diag(rs) - P[,,i]) * DD
    H[,,i] = diag(m) - M[,,i]
  }
  sigma = rep(1, n)
  Y     = array(0, dim=c(m,p))
  Q     = array(0, dim=c(m,p))
  S     = apply(M, c(1,2), mean)
  S_old = matrix(0, m, p)

  step = 0
  while(1){
    step = step + 1
    funV[step] = sum(sigma/2) +
      sum(apply(M, 3, function(x) norm(diag(m)-x-S, 'F'))/(m^2*sigma))
    relChg = norm(S - S_old, 'F')/max(1, norm(S_old, 'F'))
    S_old = S;
    if(verbose){
      cat(paste0('iter',step,
                 '\n relChg',round(relChg, 4),
                 ',\n funV=',round(funV[step],4), '\n'))
    }
    if (step > 1 && relChg < eps){
      print("converged!")
      break;
    }
    if (step > max_iter){
      print(paste('reached max iteration : ', step))
      break;
    }


    # update Q
    phi = sum(1/(2*sigma*m*m))
    tmp = matrix(0, m, m)
    for (i in 1:n){
      tmp = tmp + H[,,i] / (2*phi*sigma[i]*m*m)
    }
    USigV = svd(tmp)
    U = USigV$u; Sigma = USigV$d; V = USigV$v;
    C = tau / phi
    svp = sum(Sigma > C)
    if(svp>=2){
      Sigma = Sigma[1:svp]-C;
      S = U[, 1:svp] %*% diag(Sigma) %*% t(V[,1:svp])
    }else if(svp==1){
      Sigma = Sigma[1]-C;
      S = as.matrix(U[,1]) %*% matrix(Sigma, nrow=1) %*%  as.matrix(t(V[,1]))
    }else{
      svp=1
      S = matrix(0,m,m)
    }

    # #update sigma
    sigma = apply(H, 3, function(x) norm(x-S, 'F'))
    sigma = sigma/(m^2)
  }
  # heat(S)
  # pi = irlba(t(S), 1)$v
  # Dist = as.numeric(pi)/sum(pi)
  # Dist2 = sum(pi)/as.numeric(pi)
  # sspi = sqrt(diag(Dist2))
  # spi = sqrt(diag(Dist))
  # S = (spi %*% S %*% sspi + sspi %*% t(S) %*% spi)/2
  return(list(S=S, f=funV, sigma=sigma))
  # return(list(S=S, f = funV, sigma=sigma))
}
