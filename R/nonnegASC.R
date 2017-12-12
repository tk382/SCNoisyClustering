nonnegASC = function(B){
  n = nrow(B); m = ncol(B)
  A = t(matrix(rep(1:m, n), ncol=n))
  B_sort = t(apply(B, 1, function(x) sort(x, decreasing=TRUE)))
  cum_B = rowCumsums(B_sort)
  sigma = B_sort - (cum_B-1)/A
  tmp = apply((sigma>0), 2, as.numeric)
  idx = rowSums(tmp)
  tmp = B_sort-sigma
  sigma = diag(tmp[, idx])
  sigma = matrix(rep(sigma, m), nrow=length(sigma))
  X = pmax(B-sigma, 0);
  return(X)
}
