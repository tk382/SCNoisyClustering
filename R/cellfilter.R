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
cellfilter = function(X){
  # zeros = colSums(X==0)
  # remove = which(zeros > (p1 * ncol(X)))
  # remove = c(remove, which(zeros < (p2 * ncol(X))))
  # if(length(remove)>0){X = X[,remove]}
  # return(X)
  nUMI = colSums(X); #hist(nUMI)
  nGene = colSums(X>0);# hist(nGene)
  det.rate = nGene / nrow(X);# hist(det.rate)
  mito.genes = grep("^MT-", rownames(X), value=T)
  if(length(mito.genes)>0){
    percent.mito = colSums(X[mito.genes,]) / nUMI
  }else{
    percent.mito = rep(0, ncol(X))
  }
  keep.cells = nGene > 500 &  percent.mito < 0.1
  X = X[,keep.cells]
  nUMI = nUMI[keep.cells]
  det.rate = det.rate[keep.cells]
  nGene = nGene[keep.cells]
  percent.mito = percent.mito[keep.cells]
  det.rate = colSums(X>0) / nrow(X)
  return(list(X=X, nUMI = nUMI, nGene = nGene, percent.mito = percent.mito, det.rate = det.rate))
}
