SLSL_ref = function(X, ref, numClust){
  int = base::intersect(rownames(ref), rownames(X))
  if(length(int)<100){
    stop("Not enough genes overlap with the reference data set")
  }

  # take subset of the data and
  # the reference set with overlapping genes
  ind1 = match(int, rownames(X))
  X = X[ind1, ]
  ind2 = match(int, rownames(ref))
  ref = ref[ind2, ]

  projection = t(t(X) %*% as.matrix(ref))/nrow(X)
  projection = scale(projection^4)

  res = SLSL(projection, numClust = numClust, log=F)
  out = list(SLSL_output = res, projection = projection, numGenes=length(int))
  return(out)
}
