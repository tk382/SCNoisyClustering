#' Remove abnormal cells and return the summary of each cell including the UMI count, gene count, percent of mito-gene, and detection rate.
#'
#' @param X log-transformed expression level matrix with cells at columns, genes at rows
#' @param genenames : gene names for each row of X
#' @export
explore_data = function(X, genenames){
  nUMI = Matrix::colSums(X)
  nGene = Matrix::colSums(X>0)
  det.rate = nGene / nrow(X)
  mito.genes = grep("^MT-", genenames)
  percent.mito = rep(0, ncol(X))
  if(length(mito.genes)>0){
    percent.mito = Matrix::colSums(X[mito.genes,]) / nUMI
  }
  par(mfrow = c(1, 3))
  plot(nUMI ~ det.rate, xlab="detection rate", ylab="nUMI", main="nUMI and detection rate")
  boxplot(det.rate, main = "detection rate")
  boxplot(percent.mito, main = "proportion of mito-genes")
  summary = data.frame(nUMI = nUMI, det.rate = det.rate, percent.mito = percent.mito)
  return(summary)
}
