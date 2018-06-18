RMSC = function(P, tau, gamma, lambda, mu=1e-3, rho = 1.9, max_iter=100, eps=1e-9, verbose=FALSE){
  dims = dim(P)
  funV = rep(0,max_iter)
  m = dims[1]; p = dims[2]; n = dims[3];
  if(m!=p){
    stop('input matrix P must be square transition matrix')
  }
  E     = array(rnorm(m*p*n), dim = c(m,p,n))
  W     = array(0, dim=c(m,p,n))
  Z     = array(0, dim=c(m,p))
  Y     = array(0, dim=c(m,p))
  B     = array(0, dim=c(m,p,n))
  Q     = array(0, dim=c(m,p))
  R     = array(0, dim=c(m,p))
  S     = array(0, dim=c(m,p))
  S_old = matrix(rnorm(m*p), m, p)
  e     = rep(1, m)

  step = 0
  while(1){
    step = step + 1
    max_inf_norm = -1
    for (i in 1:n){
      Pi = P[,,i]
      Ei = E[,,i]
      diff = Pi-Ei-S
      inf_norm_of_i = norm(diff, 'I')
      max_inf_norm = max(max_inf_norm, inf_norm_of_i)
    }
    funV[step] = sum(svd(S)$d) + norm(S, 'O')+sum(lambda^2*apply(E, 3, function(x) norm(x, 'F')))
    relChg = norm(S-S_old, 'F')/max(1, norm(S_old, 'F'))
    S_old = S;
    tmp1 = S-Q;
    tmp2 = S-R;
    max_inf_norm2 = norm(tmp1, 'I')
    max_inf_norm3 = norm(tmp2, 'I')
    if(verbose){
      cat(paste0('iter',step,': \n Pi-E-S=',round(max_inf_norm,4),
            ',\n relChg=',round(relChg,4), ',\n mu=',round(mu,4),
            ', \n S-Q=', round(max_inf_norm2,4),
            ', \n S-R=', round(max_inf_norm3, 4),
            ',\n funV=',round(funV[step],4), '\n'))
    }
    if (step > 1 && max_inf_norm < eps && relChg<eps){
      break;
    }
    if (step > max_iter){
      print(paste('reached max iteration : ', step))
      break;
    }

    #update S
    B = 1/(n+2) * (Q + R -Z/mu - Y/mu + apply(E-P-W/mu, c(1,2),sum))
    S = nonnegASC(B)

    #update Q
    M = S+Y/mu
    C = tau/mu
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

    #update R
    R = pmax(S+Z/mu - gamma/mu, 0) + pmin(S+Z/mu - gamma/mu, 0)

    #update lambda
    for (i in 1:n){
      norms = apply(E, 3, function(x) norm(x, 'F'))
      inversenorms = 1/(norms^2)
      lambda = inversenorms / sum(inversenorms)
    }
    # Dmat = diag(n)
    # dvec = apply(P, 3, function(x) norm(x-S, 'F'))
    # Amat = cbind(rep(1,n), diag(n))
    # bvec = c(1, rep(0, n))
    # qp = solve.QP(Dmat, -dvec, Amat, bvec, meq=1)
    # lambda = qp$solution

    #update Ei and Wi
    for (i in 1:n){
      E[,,i] = mu/(2*lambda[i]^2 + mu) * (P[,,i]-S-W[,,i]/mu)
      W[,,i] = W[,,i] + mu * (S+E[,,i]-P[,,i])
    }

    #update Z
    Z = Z+mu*(S-R)

    #update Y
    Y = Y+mu*(S-Q)

    #update mu
    mu = 1
    # mu = min(rho*mu, 1e+10)
  }
  #ggplot(melt(S), aes(x=X1, y=X2, fill=value)) + geom_tile() +
  #  scale_color_gradient()+ggtitle(step)
  pi = irlba(t(S), 1)$v
  Dist = as.numeric(pi)/sum(pi)
  Dist2 = sum(pi)/as.numeric(pi)
  sspi = sqrt(diag(Dist2))
  spi = sqrt(diag(Dist))
  S = (spi %*% S %*% sspi + sspi %*% t(S) %*% spi)/2
  return(list(S=S, E=E, f = funV, lambda=lambda))
}
