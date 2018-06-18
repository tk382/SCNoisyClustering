RMSC_with_w = function(P, numClust, mu=1e-3, rho = 1.9, max_iter=100, eps=1e-9, verbose=TRUE){
  library(quadprog)
  dims = dim(P);
  funV = rep(0,max_iter)
  m = dims[1]; p = dims[2]; n = dims[3];
  sigma = rep(1, n)
  w     = rep(1/n, n)
  if(m!=p){
    stop('input matrix P must be square transition matrix')
  }
  Q     = array(0, dim=c(m,p,n))
  Y     = array(0, dim=c(m,p,n))
  S     = array(0, dim=c(m,p))
  S_old = matrix(rnorm(m*p), m, p)

  step = 0
  while(1){
    step = step + 1
    max_inf_norm = -1
    for (i in 1:n){
      Qi = Q[,,i]
      diff = S-Qi
      inf_norm_of_i = norm(diff, 'I')
      max_inf_norm = max(max_inf_norm, inf_norm_of_i)
    }
    funV[step] = sum(sigma/2)  + sum(w^2) +
      sum(apply(Q-P, 3, function(x) norm(x, 'F'))*w/(2*m*sigma))
    relChg = norm(S-S_old, 'F')/max(1, norm(S_old, 'F'))
    S_old = S;
    if(verbose){
      cat(paste0('iter',step,': \n Pi-E-S=',round(max_inf_norm,4),
                 ',\n relChg=',round(relChg,4), ',\n mu=',round(mu,4),
                 ',\n funV=',round(funV[step],4), '\n'))
    }
    if (step > 1 && max_inf_norm < eps){
      break;
    }
    if (step > max_iter){
      print(paste('reached max iteration : ', step))
      break;
    }

    #update Qi
    for (i in 1:n){
      tmp = (mu*S + P[,,i]*w[i]/(sigma[i]*m)+Y[,,i])
      tmp = tmp * (m * sigma[i] / (w[i] + m * sigma[i] * mu))
      s = svd(tmp)
      U = s$u; d = s$d; V = s$v
      d[(numClust+1):length(d)] = 0
      Q[,,i] = U %*% diag(d) %*% t(V)
    }

    #update S
    tmp = apply(Q - Y/mu, c(1,2), sum) / n;
    s = svd(tmp)
    U = s$u; d = s$d; V = s$v;
    d[(numClust+1):(length(d))] = 0
    S = U %*% diag(d) %*% t(V)


    #update sigma
    sigma = apply(Q-P, 3, function(x) norm(x, 'F')) * w
    sigma = sigma/sqrt(m)
    # norms = apply(Q-P, 3, function(x) norm(x, 'F'))*w
    # inversenorms = 1/(norms^2)
    # sigma = sum(inversenorms)/inversenorms

    #update w
    Dmat = diag(n)
    dvec = apply(Q-P, 3, function(x) norm(x, 'F'))/(2*m*sigma)
    Amat = cbind(rep(1,n), diag(n))
    bvec = c(1, rep(0, n))
    qp = solve.QP(Dmat, dvec, Amat, bvec, meq=1)
    w = qp$solution
    #print(w)

    #update Y
    for (i in 1:n){
      Y[,,i] = Y[,,i] + mu * (S-Q[,,i])
    }
    #update mu
    mu = min(rho*mu, 1e+10)

  }
  #ggplot(melt(S), aes(x=X1, y=X2, fill=value)) + geom_tile() +
    scale_color_gradient()+ggtitle(step)
  #plot(eigen(S)$values)
  pi = irlba(t(S), 1)$v
  Dist = as.numeric(pi)/sum(pi)
  Dist2 = sum(pi)/as.numeric(pi)
  sspi = sqrt(diag(Dist2))
  spi = sqrt(diag(Dist))
  S = (spi %*% S %*% sspi + sspi %*% t(S) %*% spi)/2
  return(list(S=S, f = funV, sigma = sigma))
}
