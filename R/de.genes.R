#' Remove genes (rows) with little information.
#' Genes that are not expressed in any of the samples, or most of the samples, and in some cases, genes that are expressed in most the samples, will be removed.
#'
#' @param X log-transformed expression level matrix with cells at columns, genes at rows
#' @param mean.thresh : numeric vector of length 2 for lower and upper bound for the mean expression level
#' @param dispersion.thresh : numeric vector of length 2 for lower and upper bound for log variance to mean ratio
#' @return expression level matrix containing only highly variable genes
#' @export
de.genes = function(X, genenames, cluster, top.n = 100, plot = 6){
  p = rep(0, nrow(X))
  for (i in 1:nrow(X)){
    p[i] = kruskal.test(X[i,], as.factor(out$result))$p.value
  }
  sorted_p = sort(-log10(p), index.return=TRUE, decreasing = TRUE)
  sorted_type = sort(cluster, index.return=TRUE)$ix
  par(mfrow = c(2,3))
  numClust = length(unique(cluster))
  for (i in 1:plot){
    plot(X[sorted_p$ix[i], sorted_type],
         col = rainbow(numClust)[as.factor(cluster[sorted_type])],
         main = genenames[sorted_p$ix[i]], cex=0.5,
         ylab = "pre-processed log expression", xlab = "cells")
  }
  return(data.frame(de_genes = genenames[sorted_p$ix[1:top.n]],
                    log10p = sorted_p$x[1:top.n]))
}
