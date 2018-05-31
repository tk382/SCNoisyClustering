#### Problem set up ####
par(mar=c(5,5,5,3))
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
library(SIMLR)
library(igraph) #for nmi

library(SC3)
library(SingleCellExperiment)
library(scater)
library(pcaReduce)
library(Rtsne)
set.seed(1)
R.utils::sourceDirectory("R/", modifiedOnly=FALSE)
Rcpp::sourceCpp('src/SCNoisyClustering.cpp')




#10X Gemcode - pbmc 3k
library(Matrix)
orig = readMM('data/pbmc3k/matrix.mtx')
dim(orig)
orig_zeros = apply(orig, 1, function(x) sum(x==0))
ind = which(orig_zeros > ((ncol(orig))-5))
orig = orig[-ind,]
orig_zeros = orig_zeros[-ind]
orig_pc1   = irlba(orig, 1)$v[,1]

orig_rs = t(t(orig) - (sum(orig_pc1^2)) *orig_pc1 %*% (t(orig_pc1) %*% t(orig)))


rs_pc1 = irlba(orig_rs, 1)$v[,1]
rs_zeros = apply(orig, 2, function(x) sum(x==0))
plot(rs_pc1 ~ rs_zeros, xlab='detection rate',
     ylab = 'pc1',
     main='count matrix after removing pc1')

plot(orig_pc1 ~ apply(orig, 2, function(x) sum(x==0)),
     xlab='detection rate', ylab='pc1',
     main='count matrix before removing pc1')


lorig = log(orig+1)
lorig_zeros = apply(lorig, 1, function(x) sum(x==0))
lorig_pc1   = irlba(lorig, 1)$v[,1]
lorig_rs = t(t(lorig) - (sum(lorig_pc1^2)) *lorig_pc1 %*% (t(lorig_pc1) %*% t(lorig)))
lrs_pc1 = irlba(lorig_rs, 1)$v[,1]
lrs_zeros = apply(lorig, 2, function(x) sum(x==0))
plot(lrs_pc1 ~ lrs_zeros, xlab='detection rate',
     ylab='pc1',
     main='log count matrix after removing pc1')

plot(lorig_pc1 ~ apply(lorig, 2, function(x) sum(x==0)),
     xlab='detection rate',
     ylab='pc1',
     main='log count matrix before removing pc1')
