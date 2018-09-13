#' Remove abnormal cells and return the summary of each cell including the UMI count, gene count, percent of mito-gene, and detection rate.
#'
#' @param X log-transformed expression level matrix with cells at columns, genes at rows
#' @param nUMI.thresh : numeric vector of length 2 for lower and upper bound for nUMI
#' @param det.rate.thresh : numeric vector of length 2 for lower and upper bound for detection rate
#' @param percent.mito.thresh : numeric vector of length 2 for lower and upper bound for the proportion of mito-genes
#'
#' @export
cell.filter = function(summary, nUMI.thresh, det.rate.thresh, percent.mito.thresh){
    ind = which(summary[,1] > nUMI.thresh[1] &
                  summary[,1] < nUMI.thresh[2] &
                  summary[,2] > det.rate.thresh[1] &
                  summary[,2] < det.rate.thresh[2] &
                  summary[,3] > percent.mito.thresh[1] &
                  summary[,3] < percent.mito.thresh[2])
    return(ind)
  }
