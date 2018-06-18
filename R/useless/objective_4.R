objective_4 = function(P, tau, gamma, r, mu=1e-3, eps = 1e-6, rho=1.9, max_iter=100, verbose=FALSE){
  m = p  = dim(P)[[1]]
  n      = dim(P)[[3]]
  S      = apply(P, c(1,2), sum)
  S_old  = matrix(0, m, p)
  Q      = matrix(0,m,p)
  D_half = 1/sqrt(rowSums(S))
  Y      = matrix(0, m,p)
  DD     = matrix(rep(D_half, nrow(S)), nrow=nrow(S))
  DD     = DD * t(DD)
  Ltilde = (S-diag(m)) * DD
  J      = irlba(Ltilde)$v[,1:4]
  w      = rep(1/n, n)
  step   = 0
  funV   = rep(0, max_iter)
  each_term = matrix(0, max_iter, 4)
  while(1){
    step = step + 1
    max_inf_norm = norm(Q-S, 'I')
    funV[step] = tau*sum(svd(S)$d) -
      2*gamma*sum(diag(t(J) %*% Ltilde %*% J)) +
      sum(w*apply(P, 3, function(x) norm(x-S, 'F'))) + sum(w^2)/2
    relChg = norm(S-S_old, 'F')/max(1, norm(S_old, 'F'))
    S_old = S;
    if (step > max_iter){
      print(paste('reached max iteration : ', step-1))
      break;
    }
    if(verbose){
      cat(paste0('iter',step,': \n S-Q=',round(max_inf_norm,4),
                 ',\n relChg=',round(relChg,4), ',\n mu=',round(mu,4),
                 ',\n funV=',round(funV[step],4), '\n'))
    }
    each_term[step,] = c(tau*sum(svd(S)$d),
                         gamma*sum(diag(t(J) %*% diag(D_half) %*% (diag(m) - S) %*% diag(D_half) %*% J)),
                         sum(w*apply(P, 3, function(x) norm(x-S, 'F'))/2),
                         sum(w^2)/2)
    if (step > 1 && max_inf_norm < eps && relChg < eps){
      print("converged!")
      break;
    }

    #update Q
    C = tau / (mu+1)
    Jtilde = D_half * J
    M = gamma * (Jtilde %*% t(Jtilde)) + Y + mu*S
    for (i in 1:n){
      M = M + w[i] * P[,,i]
    }
    M = M/(mu+1)
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

    #update J
    D_half = 1/sqrt(rowSums(S))
    DD = matrix(rep(D_half, nrow(S)), nrow=nrow(S))
    DD = DD * t(DD)
    Ltilde = (S-diag(m)) * DD
    J = svd(Ltilde)$v[,1:r]

    #update w
    Dmat = diag(n)
    dvec = apply(P, 3, function(x) norm(x-S, 'F'))^2
    dvec = dvec/2
    Amat = cbind(rep(1,n), diag(n))
    bvec = c(1, rep(0, n))
    qp = solve.QP(Dmat, -dvec, Amat, bvec, meq=1)
    w = qp$solution

    #update Y
    Y = Y+mu*(S-Q)

    #update mu
    mu = min(rho*mu, 1e+10)
  }
  ggplot(melt(S), aes(x=X1, y=X2, fill=value)) + geom_tile() +
    scale_fill_gradient()+ggtitle(step)
  pi = svd(t(S))$v[,1]
  Dist = as.numeric(pi)/sum(pi)
  Dist2 = sum(pi)/as.numeric(pi)
  sspi = sqrt(diag(Dist2))
  spi = sqrt(diag(Dist))
  S = (spi %*% S %*% sspi + sspi %*% t(S) %*% spi)/2
  return(list(S=S, f = funV, w=w, J=J, eachterm = each_term))
}
