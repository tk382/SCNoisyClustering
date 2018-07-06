genefilter = function(X, p1=0.9, p2=0){
  zeros = apply(X, 1, function(x) sum(x==0))
  remove = which(zeros > (p1 * ncol(X)))
  if(length(remove)>0){ X = X[-remove,]}
  zeros = apply(X, 1, function(x) sum(x==0))
  remove = which(zeros < (p2 * ncol(X)))
  if(length(remove)>0){X = X[-remove,]}
  return(X)
}
