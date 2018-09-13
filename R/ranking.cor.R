ranking.cor = function(X){
  X.rank = apply(X, 2, rank, ties.method="min")
  return(cor(X.rank))
}
