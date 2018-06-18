new_sparse = function(P, c, mu=1e-3, rho = 1.9, max_iter=100, eps=1e-9, verbose=TRUE){
  dims = dim(P);
  funV = rep(0,max_iter)
  m = dims[1]; p = dims[2]; n = dims[3];
  if(m!=p){
    stop('input matrix P must be square transition matrix')
  }
  Y     = array(0, dim=c(m,p))
  Q     = array(0, dim=c(m,p))
  S     = array(0, dim=c(m,p))
  S_old = matrix(rnorm(m*p), m, p)
  w     = rep(1/n,n)
  step = 0
  while(1){
    step = step + 1
    max_inf_norm = -1
    # for (i in 1:n){
    #   diff = P[,,i]-S
    #   inf_norm_of_i = norm(diff, 'I')
    #   max_inf_norm = max(max_inf_norm, inf_norm_of_i)
    # }
    max_inf_norm = norm(Q-S, 'F')
    funV[step] = norm(S, 'O')*c +
      sum(apply(P, 3, function(x) norm(x-S, 'F'))*w) +
      sum(w^2)

    relChg = norm(S-S_old, 'F')/max(1, norm(S_old, 'F'))
    S_old = S;
    if(verbose){
      cat(paste0('iter',step,': \n |S-Q|=',round(max_inf_norm,4),
                 ',\n relChg=',round(relChg,4), ',\n mu=',round(mu,4),
                 ',\n funV=',round(funV[step],4), '\n'))
    }
    if (step > 1 && max_inf_norm < eps){
      break;
    }
    if (step > max_iter){
      print(paste('reached max iteration : ', max_iter))
      break;
    }

    #update Q
    tmp = S+Y/mu
    for (i in 1:n){
      tmp = tmp + 2*P[,,i] * w[i]
    }
    tmp = tmp * (mu+2)/mu
    Q = pmax(tmp-c/(mu+2), 0) + pmin(tmp-c/(mu+2))
    Q = Q / rowSums(Q)

    #update S
    tmp = Q - Y / mu;
    s = svd(tmp)
    U = s$u; d = s$d; V = s$v;
    d[(numClust+1):(length(d))] = 0
    S = U %*% diag(d) %*% t(V)
    S = S / rowSums(S)


    #update w
    Dmat = diag(n)
    dvec = -apply(P, 3, function(x) norm(x-S, 'F'))
    Amat = cbind(rep(1,n), diag(n))
    bvec = c(1, rep(0, n))
    qp = solve.QP(Dmat, dvec, Amat, bvec, meq=1)
    w = qp$solution

    #update Y
    Y = Y + mu * (S-Q)

    #update mu
    mu = min(rho*mu, 1e+10)

  }
  ggplot(melt(S), aes(x=X1, y=X2, fill=value)) + geom_tile() +
    scale_color_gradient()+ggtitle(step)
  pi = irlba(t(S), 1)$v
  Dist = as.numeric(pi)/sum(pi)
  Dist2 = sum(pi)/as.numeric(pi)
  sspi = sqrt(diag(Dist2))
  spi = sqrt(diag(Dist))
  S = (spi %*% S %*% sspi + sspi %*% t(S) %*% spi)/2
  return(list(S=S, f = funV, w = w))
}

