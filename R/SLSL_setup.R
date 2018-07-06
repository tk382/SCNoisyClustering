#### Problem set up ####
SLSL_setup = function(){
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
  library(broom)
  library(abind)


  set.seed(1)
  R.utils::sourceDirectory("R/", modifiedOnly=FALSE)
  Rcpp::sourceCpp('src/SCNoisyClustering.cpp')
  return(NA)
}
