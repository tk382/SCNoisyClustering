adaptive_lasso = function(S, mu, eta = 1, maxiter=10){
  n = nrow(S)
  Xi = S
  Gamma = Qi = dP = list(1:n-1)
  for (ii in 1:(n-1)){
    Gamma[[ii]] = matrix(0, n-ii, n)
  }

  for (ii in 1:(n-1)){
    const_mat = cbind(matrix(0, n-ii, ii-1), matrix(1, n-ii, 1), -diag(n-ii))
    Qi[[ii]] = const_mat %*% S
    dP[[ii]] = rowSums(Qi[[ii]]^2)
  }
  DelDel = n*diag(n) - matrix(1,n,n)
  Gammai = Gamma; Q = Qi
  iter = 0; err = 10

  while ((err > 0.1) & (iter<maxiter)){
    iter = iter + 1

    #update X
    const_mat = cbind(matrix(0, n-ii, ii-1), matrix(1, n-ii, 1), -diag(n-ii))
    int_mid1 = t(const_mat) %*% Qi[[ii]]
    int_mid2 = t(const_mat) %*% Gamma[[ii]]
    for(ii in 2:(n-1)){
      const_mat = cbind(matrix(0, n-ii, ii-1), matrix(1, n-ii, 1), -diag(n-ii))
      int_mid1 = int_mid1 + t(const_mat) %*% Qi[[ii]]
      int_mid2 = int_mid2 + t(const_mat) %*% Gamma[[ii]]
    }
    X = solve(eta * DelDel + 2*diag(n)) %*% (2*S + eta*int_mid1 - int_mid2)
    #heat(X)
    #X = nonnegASC(X)

    #update Q
    err_Q = 0
    for (ii in 1:(n-1)){
      const_mat = cbind(matrix(0, n-ii, ii-1), matrix(1, n-ii, 1), -diag(n-ii))
      tjk = Gamma[[ii]] + eta*const_mat %*% X
      tjk = tjk/mu
      tjknorm = sqrt(rowSums(abs(tjk)^2))
      Q[[ii]]= as.numeric(tjknorm>(1/dP[[ii]]))*
        ((mu*(tjknorm-(1/dP[[ii]]))/(eta*pmax(tjknorm,1e-20)))*tjk);
      err_Q = err_Q + norm(Q[[ii]] - Qi[[ii]], 'F')^2
    }

    #update Gamma

    err_Gamma = 0
    for (ii in 1:(n-1)){
      const_mat = cbind(matrix(0, n-ii, ii-1), matrix(1, n-ii, 1), -diag(n-ii))
      Gamma[[ii]] = Gammai[[ii]] + 0.1*eta*const_mat%*%X-Q[[ii]]
      err_Gamma = err_Gamma + norm(Gamma[[ii]] - Gamma[[ii]], 'F')^2
    }
    err = norm(Xi-X, 'F')^2 + err_Gamma
    Xi = X; Qi = Q; Gammai = Gamma
  }
  X = max(X) - X
  X = (X + t(X))/2
  return(X)
}
