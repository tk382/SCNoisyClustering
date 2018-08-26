#' Estimate the optimal number of clusters using eigengap
#'
#' @param S estimated similarity matrix.
#' @examples
#' #create positive definite symmetric matrix
#' X = matrix(rnorm(50), nrow = 10)
#' S= t(X) %*% X
#' getClustNum(S)
#' @export
getClustNum = function(S){
  s = svd(S)
  eigengap = abs(diff(s$d))
  c = max(3, which.max(eigengap[-1])+1)
  return(c)
}
