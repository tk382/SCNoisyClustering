#' Returns scatterplot of detection rate and the first pc of the log-transformed expression level matrix
#' and the matrix with the detection rate regressed out
#' Users can decide whether to use the raw matrix or the residual matrix based on the visualization
#'
#' @param X log-transformed expression level matrix with cells at columns, genes at rows
#' @param det.rate detection rate : usually colSums(X>0) / nrow(X) where X is count matrix.
#' @examples
#'
#' X = SCNoisyClustering::yan
#' det.rate = colSums(X>0) / nrow(X)
#' out = correct_detection_rate(X, det.rate)
#' out$plot #check if there's a distinct relationship
#' X = out$residual #if you decide to regress out the detection rate
#'
#' @export
correct_detection_rate = function(X, det.rate){
  if(sum(X<0) ==0){
    warning("X must be a log-transformed matrix")
  }
  pc1 = irlba(X,1)$v[,1]
  det.rate2 = qr(cbind(rep(1,length(det.rate)), det.rate))

  residual = t(qr.resid(det.rate2, t(logX)))

  newpc1 = irlba(residual,1)$v[,1]
  par(mfrow = c(1,2))
  plt = plot(pc1 ~ det.rate, xlab="detection rate", ylab="pc1", main="before correction")
  plt2 = plot(newpc1 ~ det.rate, xlab="detection rate", ylab="pc1", main="after correction")
  return(list(plt = plt, plt2=plt2, residual = residual))
}
