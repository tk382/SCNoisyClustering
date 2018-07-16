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
lorig = log(orig + 1)
lorig_pc1 = irlba(lorig, 1)$v[,1]

plot(lorig_pc1 ~ log(orig_zeros))
plot(lorig_pc1 ~ orig_zeros)

x = orig_zeros - mean(orig_zeros)
lorig_rs = t(t(scale(lorig)) - x %*% (t(x) %*% t(scale(lorig))/(sum(x^2))))


rs_pc1 = irlba(lorig_rs, 1)$v[,1]
rs_zeros = apply(orig, 2, function(x) sum(x==0))
plot(rs_pc1 ~ rs_zeros, xlab='detection rate',
     ylab = 'pc1',
     main='log count matrix after removing pc1')

plot(lorig_pc1 ~ apply(lorig, 2, function(x) sum(x==0)),
     xlab='detection rate', ylab='pc1',
     main='log count matrix before removing pc1')


Rtsne(lorig)
