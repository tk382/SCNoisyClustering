genefilter = function(X){
  zeros = apply(X, 1, function(x) sum(x==0))
  remove = which(zeros > (0.95 * ncol(X)))
  if(length(remove)>0){ X = X[-remove,]}
  return(X)
}
