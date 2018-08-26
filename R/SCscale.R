#' Scale each gene to have mean 0 and variance 1
#'
#'
#' @param X log-transformed expression level matrix
#'
#'
#' @export
SCscale = function(X){
  if(sum(X<0)==0){
    warning("log-transform X before scaling")
  }
  t(scale(t(X)))
}
