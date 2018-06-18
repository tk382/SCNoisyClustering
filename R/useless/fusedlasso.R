fusedlasso = function(Y, k){
  wts = rep(1, nrow(S)-1)
  phi = rep(1,nrow(S))
  grp.reg = doGFLars(Y, k, phi, wts)
  p = ncol(S); n = nrow(S)
  bkp <- grp.reg$bkp
  order.bkp <- order(bkp)
  bkp <- bkp[order.bkp]
  change.p <- bkp + 1

  delta2 <- matrix(0, p - 1, n)
  delta2[bkp, ] <- grp.reg$value[[k]][order.bkp, ]

  temp <- cumsum(phi[p : 1]^2)[p : 1][-1]
  delta2.up <- delta2 * wts
  delta1 <- t(phi) %*% Y / p - t(temp) %*% delta2.up / p

  delta <- rbind(delta1, delta2.up)
  theta <- delta
  for (i in 2 : p) theta[i, ] <- theta[i - 1, ] + delta[i, ]
  return((theta + t(theta))/2)
}
