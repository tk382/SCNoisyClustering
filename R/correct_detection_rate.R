#' Returns scatterplot of detection rate and the first pc of the log-transformed expression level matrix
#' and the matrix with the detection rate regressed out
#' Users can decide whether to use the raw matrix or the residual matrix based on the visualization
#'
#' @param X log-transformed expression level matrix with cells at columns, genes at rows
#' @examples
#'
#' X = SCNoisyClustering::yan
#' det.rate = colSums(X>0) / nrow(X)
#' out = correct_detection_rate(X, det.rate)
#' out$plot #check if there's a distinct relationship
#' X = out$residual #if you decide to regress out the detection rate
#'
#'
correct_detection_rate = function(X, det.rate){
  pc1 = irlba(X,1)$v[,1]
  det.rate2 = qr(cbind(rep(1,length(det.rate)), det.rate))
  plt = plot(pc1 ~ det.rate)
  residual = t(qr.resid(det.rate2, t(logX)))
  return(list(plt = plt, residual = residual))
}
