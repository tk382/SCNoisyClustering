nuclear_objective = function(P, tau, mu=1e-3, rho = 1.9, max_iter=100, eps=1e-9, verbose=FALSE){
  dims = dim(P)
  funV = rep(0,max_iter)
  m = dims[1]; p = dims[2]; n = dims[3];
  if(m!=p){
    stop('input matrix P must be square transition matrix')
  }
  w     = rep(1/n, n)
  Y     = array(0, dim=c(m,p))
  Q     = array(0, dim=c(m,p))
  S     = apply(P, c(1,2), mean)
  S_old = matrix(rnorm(m*p), m, p)

  step = 0
  while(1){
    step = step + 1
    max_inf_norm = norm(S-Q, 'I')
    funV[step] = tau * sum(svd(S)$d) +
      sum(w*apply(P, 3, function(x) norm(x-Q, 'F'))) + sum(w^2)
    relChg = norm(S - S_old, 'F')/max(1, norm(S_old, 'F'))
    S_old = S;
    if(verbose){
      cat(paste0('iter',step,': \n S-Q=',round(max_inf_norm,4),
                 ',\n relChg=',round(relChg,4), ',\n mu=',round(mu,4),
                 ',\n funV=',round(funV[step],4), '\n'))
    }
    if (step > 1 && max_inf_norm<eps && relChg<eps){
      print("converged!")
      break;
    }
    if (step > max_iter){
      print(paste('reached max iteration : ', step))
      break;
    }


    # update Q
    M = mu * S + Y
    for (i in 1:n){
      M = M + w[i]*P[,,i]
    }
    M = M/(mu+2)
    C = tau/(mu+2)
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

    #update S
    S = nonnegASC(Q-Y/mu)



    # update w
    Dmat = diag(n)
    dvec = apply(P, 3, function(x) norm(x-S, 'F'))
    dvec = dvec^2
    Amat = cbind(rep(1,n), diag(n))
    bvec = c(1, rep(0, n))
    qp = solve.QP(Dmat, -dvec, Amat, bvec, meq=1)
    w = qp$solution

    #update Y
    Y = Y+mu*(S-Q)

    #update mu
    mu = min(rho*mu, 1e+10)
    # mu = 1
  }
  #ggplot(melt(S), aes(x=X1, y=X2, fill=value)) + geom_tile() +
  #  scale_color_gradient()+ggtitle(step)
  pi = irlba(t(S), 1)$v
  Dist = as.numeric(pi)/sum(pi)
  Dist2 = sum(pi)/as.numeric(pi)
  sspi = sqrt(diag(Dist2))
  spi = sqrt(diag(Dist))
  S = (spi %*% S %*% sspi + sspi %*% t(S) %*% spi)/2
  return(list(S=S, f = funV, w=w))
}
