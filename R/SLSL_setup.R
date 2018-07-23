#### Problem set up ####
SLSL_setup = function(){
  library(data.table)
  library(inline)
  library(matrixStats)
  library(quadprog)
  library(irlba)
  library(ggplot2)
  library(dplyr)
  library(reshape)
  library(caret)
  library(fossil)
  library(pracma)
  library(igraph)
  library(Rtsne)
  library(gplots)
  library(broom)
  library(abind)
  library(stargazer)
  library(scatterplot3d)
  library(diceR)
  library(parallel)


  set.seed(1)
  R.utils::sourceDirectory("R/", modifiedOnly=FALSE)
  Rcpp::sourceCpp('src/SCNoisyClustering.cpp')
  return(NA)
}
