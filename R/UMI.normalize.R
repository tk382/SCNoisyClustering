#' Scale each gene to have mean 0 and variance 1
#'
#'
#' @param X Raw count matrix
#'
#'
#' @export
UMI_normalize = function(X){
  nUMI = Matrix::colSums(X)
  med = median(nUMI)
  X = t(t(X) / nUMI) * med
  return(X)
}
