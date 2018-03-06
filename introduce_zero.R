##introduce 0's in a simulated data

#### Problem set up ####
library(matrixStats)
library(irlba)
library(ggplot2)
library(KRLS)
library(dplyr)
library(reshape)
library(caret)
library(fossil)
library(gridExtra)
library(pracma)
library(igraph)
library(parallel)
library(scatterplot3d)
library(microbenchmark)
library(inline)
library(quadprog)
library(SIMLR)
source("../SIMLR/R/tsne.R")
set.seed(1)
R.utils::sourceDirectory("R/", modifiedOnly=FALSE)
Rcpp::sourceCpp('src/SCNoisyClustering.cpp')

load('../data/Test_3_Pollen.RData')
X        = Test_3_Pollen$in_X
truth    = Test_3_Pollen$true_labs
numClust = Test_3_Pollen$n_clust
ind      = sort(truth$V1, index.return=TRUE)$ix
truth    = truth$V1[ind]
rm(Test_3_Pollen)
X = X[, ind]

zeros = seq(10000,50000,by=10000)
out = data.frame(ri=rep(0, length(zeros)), ari=rep(0, length(zeros)), nmi=rep(0, length(zeros)))
for (i in 1:length(zeros)){
  numzeros             = zeros[i]
  rowind               = sample(1:nrow(X), size=numzeros, replace = TRUE)
  colind               = sample(1:ncol(X), size=numzeros, replace = TRUE)
  newX                 = X;
  newX[rowind, colind] = 0
  allk                 = seq(10, 30, by=5);
  sigma                = seq(0.2,1, by=0.2);
  P                    = multiple_kernel_tae(t(X), allk, sigma)
  n                    = nuclear_objective_c(P, 1)
  clust_n_tsne         = tsne_spectral(n$S, numClust)
  ri                   = rand.index(truth, clust_n_tsne)
  ari                  = adj.rand.index(truth, clust_n_tsne)
  nmi                  = compare(truth, clust_n_tsne, method='nmi')
  out[i, ]             = c(ri, ari, nmi)
}
out
