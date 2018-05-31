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

#### set up data ####
load('data/1_Kolod.RData')
X1        = Test_2_Kolod$in_X;
truth1    = Test_2_Kolod$true_labs
numClust1 = Test_2_Kolod$n_clust;
rm(Test_2_Kolod)

P = corr_kernel_c(t(X1), seq(10,30,by=5), seq(2,1,by=-0.1), 23)
S = sparse_scaledlasso_c(P, 5, 0, verbose=TRUE)
getClustNum(S$S) #correct

rm(X1, truth1, numClust1, P, S); gc()

load('data/2_Pollen.RData')
X2        = Test_3_Pollen$in_X
truth2    = Test_3_Pollen$true_labs
numClust2 = Test_3_Pollen$n_clust
rm(Test_3_Pollen)

X2 = genefilter(X2)
P = corr_kernel_c(t(X2), seq(10,30,by=5), seq(2,1,by=-0.1), ncol(X2)/20)
S = sparse_scaledlasso_c(P, 5, 0, verbose=TRUE)
getClustNum(S$S) #incorrect : 3 instead of 11

rm(X2, truth2, numClust2, P, S); gc()

load('data/3_Usoskin.RData')
X3        = Test_4_Usoskin$in_X
truth3    = Test_4_Usoskin$true_labs
numClust3 = Test_4_Usoskin$n_clust
rm(Test_4_Usoskin)

X3 = genefilter(X3)
P = corr_kernel_c(t(X3), seq(10,30,by=5), seq(2,1,by=-0.1), ncol(X3)/20)
S = sparse_scaledlasso_c(P, 5, 1, verbose=TRUE)
getClustNum(S$S) #incorrect : 15 instead of 11
res = tsne_spectral(S$S, 4)
adj.rand.index(res$cluster, truth3$V1)

Usostsne = Rtsne(S$S, dims=3, perplexity=30, verbose=TRUE, max_iter=500)
Y = Usostsne$Y
plot(Y[,1], Y[,2], col = rainbow(4)[res$cluster])
scatterplot3d(Y[,1], Y[,2], Y[,3], color=rainbow(4)[res$cluster])

rm(X3, S, res, Usostsne, Y)
gc()

load('data/4_Buettener.RData')
X4        = Test_1_mECS$in_X
truth4    = Test_1_mECS$true_labs; truth4 = data.frame(V1=truth4)
numClust4 = Test_1_mECS$n_clust
rm(Test_1_mECS)

X4 = genefilter(X4)
P = corr_kernel_c(t(X4), seq(10,30,by=5), seq(2,1,by=-0.1), ncol(X4)/20)
S = sparse_scaledlasso_c(P, 5, 0, verbose=TRUE)
getClustNum(S$S) #correct
res = tsne_spectral(S$S, 3)
adj.rand.index(res$cluster, truth4$V1)

rm(X4, P , S, res)


load('data/5_Yan.rda')
X5        = as.matrix(yan)
truth5    = as.numeric(ann$cell_type1); truth5 = data.frame(V1=truth5)
numClust5 = length(unique(truth5$V1))
rm(ann, yan)

X5 = log(X5 + 1)
X5 = genefilter(X5)
P = corr_kernel_c(t(X5), seq(10,30,by=5), seq(2,1,by=-0.1), 10)
S = sparse_scaledlasso_c(P, 5, 0, verbose=TRUE)
getClustNum(S$S) #got 16..
res = tsne_spectral(S$S, numClust5)
adj.rand.index(res$cluster, truth5$V1)
rm(X5, P, S, res)

load('data/6_Treutlein.rda')
X6        = as.matrix(treutlein)
truth6    = as.numeric(colnames(treutlein))
ind       = sort(truth6, index.return=TRUE)$ix
X6        = X6[,ind]
truth6    = truth6[ind]; truth6 = data.frame(V1=truth6)
numClust6 = length(unique(truth6$V1))
rm(treutlein)

X6 = log(X6 + 1)
X6 = genefilter(X6)
P = corr_kernel_c(t(X6), seq(10,30,by=5), seq(2,1,by=-0.1), 10)
S = sparse_scaledlasso_c(P, 5, 0, verbose=TRUE)
getClustNum(S$S) #incorrect got 3 instead of 5
res = tsne_spectral(S$S, numClust6)
adj.rand.index(res$cluster, truth6$V1)
rm(X6, P, S, res)
gc()



X7        = read.table('data/7_GSE75748_sc_cell_type_ec_cleaned.txt', header=TRUE, stringsAsFactors = FALSE)
X7        = as.matrix(X7)
truth7    = data.frame(V1 = as.factor(sapply(strsplit(colnames(X7), '_'), function(x) x[1])))
truth7$V1 = as.integer(truth7$V1)
numClust7 = 7
X7 = genefilter(X7)
P = corr_kernel_c(t(X7), seq(10,30,by=5), seq(2,1,by=-0.1), ncol(X7)/20)
S = sparse_scaledlasso_c(P, 5, 0, verbose=TRUE)
S = network.diffusion(S$S, ncol(X7)/20, 0.7)
getClustNum(S) #got 7 after network diffusion
res = tsne_spectral(S, 7)
adj.rand.index(res$cluster, truth7$V1) #perfect clustering
rm(X7, P, S, res)




X8        = read.table('data/8_cleaned_celltype_exp2.txt', header=TRUE, stringsAsFactors = FALSE)
X8        = as.matrix(X8)
truth8    = data.frame(V1 = as.factor(sapply(strsplit(colnames(X8), '_'), function(x) x[1])))
truth8$V1 = as.integer(truth8$V1)
numClust8 = length(unique(truth8$V1))
X8 = genefilter(X8)
P = corr_kernel_c(t(X8), seq(10,30,by=5), seq(2,1,by=-0.1), ncol(X8)/20)
S = sparse_scaledlasso_c(P, 5, 0, verbose=TRUE)
S = network.diffusion(S$S, ncol(X8)/20, 0.7)
getClustNum(S) #got 7 after network diffusion
res = tsne_spectral(S, 7)
adj.rand.index(res$cluster, truth8$V1) #perfect clustering
rm(X8, P, S, res)



load('data/9_Macosko.Rdata')
X9 = Macosko$X
rownames(X9) = X9$gene
X9 = X9[-1]
X9 = as.matrix(X9)
zeros = apply(X9, 1, function(x) sum(x==0))
ind = which(zeros>(4800*0.95))
X9 = X9[-ind, ]
#truth9 = data.frame(V1=Macosko$label)
#numClust9 = length(unique(truth9$V1))
numClust9 = 3
rm(Macosko); gc()
