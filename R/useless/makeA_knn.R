makeA_knn = function(S, k){
  sorted = t(apply(S, 2, function(x) sort(x, index.return=TRUE, decreasing=TRUE)$ix))
  indmat = sorted[, 2:(k+1)]
  S2 = matrix(0, nrow(S), ncol(S))
  for (i in 1:nrow(S2)){
    S2[i,indmat[i,]] = 1
  }
  return(S2)
}
