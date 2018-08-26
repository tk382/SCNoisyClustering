#' Remove abnormal cells and return the summary of each cell including the UMI count, gene count, percent of mito-gene, and detection rate.
#'
#' @param X log-transformed expression level matrix with cells at columns, genes at rows
#' @param genenames : gene names for each row of X
#' @param minGene : lower bound of detected genes of a cell
#' @param maxGene : upper bound of detected genes of a cell
#' @param maxMitoProp : maximum possible value of the proportion of mito-gene reads per cell.
#'
#' @export
cellFilter = function(X, genenames, minGene = -Inf, maxGene = Inf, maxMitoProp = 0.1){
  nUMI = Matrix::colSums(X);
  nGene = Matrix::colSums(X>0);
  det.rate = nGene / nrow(X);
  mito.genes = grep("^MT-", genenames)
  if(length(mito.genes)>0){
    percent.mito = Matrix::colSums(X[mito.genes,]) / nUMI
  }else{
    percent.mito = rep(0, ncol(X))
  }
  keep.cells = nGene > minGene & nGene < maxGene & percent.mito < 0.1
  X = X[,keep.cells]
  nUMI = nUMI[keep.cells]
  det.rate = det.rate[keep.cells]
  nGene = nGene[keep.cells]
  percent.mito = percent.mito[keep.cells]
  det.rate = Matrix::colSums(X>0) / nrow(X)
  return(list(X=X, nUMI = nUMI, nGene = nGene, percent.mito = percent.mito, det.rate = det.rate))
}
