#' Remove genes (rows) with little information.
#' Genes that are not expressed in any of the samples, or most of the samples, and in some cases, genes that are expressed in most the samples, will be removed.
#'
#'
#' @param X data matrix with zeros
#' @param p1 upper threshold of the percentage of zeros in each row
#' @param p2 lower threshold of the percentage of zeros in each row
#' @return The data matrix with rows with too many or too few zeros removed
#' @examples
#' X = matrix(sample(0:10, size = 10000, replace = TRUE, prob = c(0.9, rep(0.1/10, 10))), nrow = 200) #create expression level matrix
#' newX = genefilter(X, 0.9, 0)
#' dim(newX)
genefilter = function(X, p1 = 0.9, p2 = 0){
  zeros = apply(X, 1, function(x) sum(x==0))
  remove = which(zeros > (p1 * ncol(X)))
  if(length(remove)>0){ X = X[-remove,]}
  zeros = apply(X, 1, function(x) sum(x==0))
  remove = which(zeros < (p2 * ncol(X)))
  if(length(remove)>0){X = X[-remove,]}
  return(X)
}
