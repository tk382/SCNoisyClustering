objective_6 = function(P, lambda, r, eps = 1e-6, max_iter=100, verbose=FALSE){
  m = p  = dim(P)[1]
  S      = nonnegASC(P)
  L      = diag(rowSums(S))-(S+t(S))/2
  S_old  = matrix(0, m, p)
  J      = svd(L)$v[,(ncol(L)-r+1) : ncol(L)]
  funV   = rep(0, max_iter)
  step   = 0

  while(1){
    step = step + 1
    funV[step]= 2*lambda*sum(diag(t(J) %*% L %*% J)) + norm(S-P, 'F')^2
    relChg = norm(S-S_old, 'F')/max(1, norm(S_old, 'F'))
    S_old  = S;
    if (step > max_iter){
      print(paste('reached max iteration : ', step-1))
      break;
    }
    if(verbose){
      cat(paste0('iter',step,
                 ',\n relChg=',round(relChg,4),
                 ',\n funV=',round(funV[step],4), '\n'))
    }
    if (step > 3 && relChg < eps){
      print("converged!")
      break;
    }

    #update S
    # V = J %*% t(J); eta = rep(0, m)
    # for (i in 1:m){
    #   temporaryvalue = A[i, ]-lambda*V[i,]/2
    #   eta[i] = get_eta(A[i,] - lambda*V[i,]/2)$minimum
    #   S[i,] = pmax(A[i,] - lambda * V[i,]/2 + eta[i], 0)
    # }
    V = J %*% t(J)
    S = nonnegASC(P+lambda*V/2)

    #update J
    J = svd(diag(S)-(t(S)+S)/2)$v[, (ncol(S)-r+1) : ncol(S)]
    heat(S)
  }
  plot(svd(S)$d)
  heat(S)
  pi = svd(t(S))$v[,1]
  Dist = as.numeric(pi)/sum(pi)
  Dist2 = sum(pi)/as.numeric(pi)
  sspi = sqrt(diag(Dist2))
  spi = sqrt(diag(Dist))
  S = (spi %*% S %*% sspi + sspi %*% t(S) %*% spi)/2
  return(list(S=S, f = funV, w=w, J=J, eachterm = each_term))
}

pmaxing = function(eta){
  abs(sum(pmax(temporaryvalue+eta, 0))-1)
}
get_eta = function(x){
  temporaryvalue = x
  optimize(pmaxing, c(-200, 0))
}

