#' Remove genes (rows) with little information.
#' Genes that are not expressed in any of the samples, or most of the samples, and in some cases, genes that are expressed in most the samples, will be removed.
#'
#' @param X log-transformed expression level matrix with cells at columns, genes at rows
#' @param mean.thresh : numeric vector of length 2 for lower and upper bound for the mean expression level
#' @param dispersion.thresh : numeric vector of length 2 for lower and upper bound for log variance to mean ratio
#' @return expression level matrix containing only highly variable genes
#' @export
gene_filter = function(X, genenames, dispersion, mean.thresh, dispersion.thresh){
  ind = which(dispersion$genemeans > mean.thresh[1] &
          dispersion$genemeans < mean.thresh[2] &
          dispersion$vmr > dispersion.thresh[1] &
          dispersion$vmr < dispersion.thresh[2])
  return(ind)
}

