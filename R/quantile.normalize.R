#' quantile normalize the gene expression level matrix
#'
#' @param X matrix to quantile-normalize
#'
#' @export
quantile.normalize = function(X){
  X.rank = apply(X, 2, rank, ties.method="min")
  X.sorted = apply(X, 2, sort)
  X.mean = apply(X.sorted, 1, mean)
  X.final = apply(X.rank, 2, function(x) X.mean[x])
  return(log(X.final+1))
}
