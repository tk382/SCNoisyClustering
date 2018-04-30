baseline_spectral = function(S, numClust){
  D = diag(apply(S,MARGIN=2,FUN=sum))
  L = diag(1/sqrt(diag(D))) %*% (D - S) %*% diag(1/sqrt(diag(D)))
  evL = eigen(L, symmetric=TRUE)
  Z   = evL$vectors[,(ncol(evL$vectors)-numClust+1):ncol(evL$vectors)]
  C = kmeans(Z, numClust, nstart = 200, iter.max=20)
  return(C$cluster)
}



