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
##MUST LOAD BEFORE COMPILING RCPP##
library(inline)
library(quadprog)

set.seed(1)
R.utils::sourceDirectory("R/", modifiedOnly=FALSE)
Rcpp::sourceCpp('src/SCNoisyClustering.cpp')


data_index = 1; data_name = 'mECS'
load('../data/Test_1_mECS.RData')
dat = Test_1_mECS
rm(Test_1_mECS); gc()
X = dat$in_X
truth = dat$true_labs
numClust = dat$n_clust

allk = seq(10,30,by=2); sigma=seq(0.2,1,by=0.2)
P = multiple_kernel_tae(t(X), allk, sigma)

res = rankinfo(P, 1, 1, numClust, mu=1, eps=1e-9, verbose=TRUE)
plot(svd(res$S)$d)
ggplot(melt(res$S), aes(x=X1, y=X2, fill=value))+geom_tile()+scale_fill_gradient()
