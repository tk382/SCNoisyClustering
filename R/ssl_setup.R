#### Problem set up ####
ssl_setup = function(){
  library(matrixStats)
  library(inline)
  library(quadprog)
  library(irlba)
  library(ggplot2)
  library(dplyr)
  library(reshape)
  library(caret)
  library(fossil)
  library(pracma)
  library(scatterplot3d)
  library(igraph) #for nmi
  library(Rtsne)
  library(diceR)
  library(gplots)

  set.seed(1)
  R.utils::sourceDirectory("R/", modifiedOnly=FALSE)
  Rcpp::sourceCpp('src/SCNoisyClustering.cpp')
  return(NA)
}
