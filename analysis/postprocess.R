#### Problem set up ####
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

library(SC3)
library(SingleCellExperiment)
library(scater)
library(CIDR)
library(pcaReduce)
library(Rtsne)


set.seed(1)
R.utils::sourceDirectory("R/", modifiedOnly=FALSE)
Rcpp::sourceCpp('src/SCNoisyClustering.cpp')



load('data/4_Buettener.RData')
X        = Test_1_mECS$in_X
truth    = Test_1_mECS$true_labs$V1
numClust = Test_1_mECS$n_clust
rm(Test_1_mECS)


P = mystery_kernel(t(X), k=10, allk_input = c(10,20), sigma_input = c(1,2))
res = sparse_scaledlasso(P, 5, 3, verbose=TRUE)
S = res$S
heat(S)
h = hclust(as.dist(max(S)-S), 'ave')
heat(S[h$order, h$order])
newS = fusedlasso(S[h$order, h$order], 10)
heat(newS)
new_orig_S = fusedlasso(S, 10)
heat(new_orig_S[h$order, h$order])
