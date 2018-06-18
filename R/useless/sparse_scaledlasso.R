sparse_scaledlasso = function(P, tau, gamma, k, mu=1e-3, rho = 1.9, max_iter=100, eps=1e-9, verbose=FALSE){
  dims = dim(P)
  funV = rep(0,max_iter)
  m = dims[1]; p = dims[2]; n = dims[3];
  if(m!=p){
    stop('input matrix P must be square transition matrix')
  }
  sigma    = rep(1, n)
  Y        = array(0, dim=c(m,p))
  Q        = array(0, dim=c(m,p))
  S        = apply(P, c(1,2), mean)
  S_old    = matrix(0, m, p)
  L1_norm  = F_norm = rep(0, max_iter)
  fro_norm = rep(0, max_iter)
  noise    = rep(0, max_iter)
  step = 0
  while(1){
    if (step > max_iter){
      print(paste('reached max iteration : ', max_iter))
      break;
    }
    step = step + 1
    max_inf_norm = norm(S-Q, 'I')
    L1_norm[step] = norm(S, 'O')
    fro_norm[step] = sum(apply(P, 3, function(x) norm(x-S, 'F'))/(m^2*sigma))
    F_norm[step] = norm(S, 'F')
    noise[step]     = sum(sigma)/2
    funV[step] = noise[step] +
                  fro_norm[step] +
                  tau * L1_norm[step] +
                  gamma * F_norm[step]

    relChg = norm(S - S_old, 'F')/max(1, norm(S_old, 'F'))
    S_old = S;
    if(verbose){
      cat(paste0('iter',step,': \n S-Q=',round(max_inf_norm,4),
                 ',\n relChg=',round(relChg,4), ',\n mu=',round(mu,4),
                 ',\n funV=',round(funV[step],4), '\n'))
    }
    if (step > 1 && max_inf_norm < eps && relChg < eps){
      print("converged!")
      break;
    }



    # update Q
    M = (mu*S + Y)
    phi = sum(1/(m^2 * sigma)) + mu + gamma
    for (i in 1:n){
      M = M + P[,,i] / (m^2 * sigma[i])
    }
    M = M/phi
    C = tau/phi
    Q = pmax(M-C, 0) + pmin(M-C, 0)

    #update S
    S = network.diffusion(S, k)
    S = nonnegASC(Q-Y/mu)




    # #update sigma
    sigma = apply(P, 3, function(x) norm(x-Q, 'F'))
    sigma = sigma/(m^2)

    #update Y
    Y = Y+mu*(S-Q)

    #update mu
    mu = min(rho*mu, 1e+10)
  }
  pi = irlba(t(S), 1)$v
  Dist = as.numeric(pi)/sum(pi)
  Dist2 = sum(pi)/as.numeric(pi)
  sspi = sqrt(diag(Dist2))
  spi = sqrt(diag(Dist))
  S = (spi %*% S %*% sspi + sspi %*% t(S) %*% spi)/2
  return(list(S=S,
              f = funV,
              sigma=sigma,
              L1_norm = L1_norm,
              F_norm = F_norm,
              fro_norm = fro_norm,
              noise = noise))
}
