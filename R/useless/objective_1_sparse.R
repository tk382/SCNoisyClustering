objective_1_sparse = function(X, numClust, allk = seq(2,20,by=2), sigma=seq(2,1,by=-0.25), P = NA, mu=NA, rho = 1.9, max_iter=100, eps=1e-9, verbose=TRUE){
  if(is.na(mu)){mu = 1e-3}
  if(length(dim(P))==0){
    if(verbose){
      print('computing kernel..')
    }
    K = multiple_kernel_new(t(X), allk, sigma, cores.ratio=1)
    l = length(allk) * length(sigma)
    P = array(0,dim=c(ncol(X),ncol(X),l))
    for (i in 1:l){
      P[,,i] = as.matrix(K[[i]])
    }
    rm(K)
  }
  #start algorithm
  library(quadprog)
  dims = dim(P);
  funV = rep(0,max_iter)
  m = dims[1]
  p = dims[2]
  n = dims[3]
  w = rep(1/n, n)
  c = rep(1, n)
  if(m!=p){
    stop('input matrix P must be square transition matrix')
  }
  E     = array(rnorm(m*p*n), dim=c(m,p,n))
  Y     = array(0, dim=c(m,p,n))
  S     = array(0, dim=c(m,p))
  S_old = matrix(rnorm(m*p), m, p)

  step = 0
  while(1){
    step = step + 1
    if (step > max_iter){
      print(paste('reached max iteration : ', max_iter))
      break;
    }
    max_inf_norm = -1
    for (i in 1:n){
      diff = P[,,i]-E[,,i]-S
      inf_norm_of_i = norm(diff, 'I')
      max_inf_norm = max(max_inf_norm, inf_norm_of_i)
    }
    funV[step] =  sum(w^2) + sum(c^2) +
      sum(apply(E, 3, function(x) norm(x, 'F'))*w)+
      sum(apply(E, 3, function(x) norm(x, 'O'))*c)
    relChg = norm(S-S_old, 'F')/max(1, norm(S_old, 'F'))
    S_old = S;
    if(verbose){
      cat(paste0('iter',step,': \n |S+Ei-Pi|=',round(max_inf_norm,4),
                 ',\n relChg=',round(relChg,4), ',\n mu=',round(mu,4),
                 ',\n funV=',round(funV[step],4), '\n'))
    }
    if (step > 1 && max_inf_norm < eps && relChg < eps){
      break;
    }

    #update S
    tmp = apply(E-P-Y/mu, c(1,2), sum) / n;
    s = svd(tmp)
    U = s$u; d = s$d; V = s$v;
    d[(numClust+1):(length(d))] = 0
    S = U %*% diag(d) %*% t(V)
    S = S / rowSums(S)

    #update Ei
    for (i in 1:n){
      tmp = P[,,i]-S+Y[,,i]/mu
      tmp = tmp / (2*w[i] + mu)
      E[,,i] = pmax(tmp - c[i]/(2*w[i]+mu), 0) + pmin(tmp - c[i]/(2*w[i]+mu), 0)
    }

    #update w
    Dmat = diag(n)
    dvec = -apply(E, 3, function(x) norm(x, 'F'))^2
    Amat = cbind(rep(1,n), diag(n))
    bvec = c(1, rep(0, n))
    qp = solve.QP(Dmat, dvec, Amat, bvec, meq=1)
    w = qp$solution


    #update c
    Dmat = diag(n)
    dvec = -apply(E, 3, function(x) norm(x, 'O'))
    Amat = cbind(rep(1,n), diag(n))
    bvec = c(n, rep(0, n))
    qp = solve.QP(Dmat, dvec, Amat, bvec, meq=1)
    c = qp$solution

    #update Y
    for (i in 1:n){
      Y[,,i] = Y[,,i] + mu * (S+E[,,i]-P[,,i])
      #ggplot(melt(Y[,,1]), aes(x=X1, y=X2, fill=value)) + geom_tile() +
      #  scale_color_gradient()+ggtitle(step)
    }

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
  return(list(S=S, f=funV[1:step], w = w))
}
