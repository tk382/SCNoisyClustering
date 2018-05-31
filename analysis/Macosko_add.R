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


final = read.table('analysis/final_benchmark_comparison.txt', header=TRUE)
newdf = data.frame(method = c('kmeans','pcareduce','SC3','SIMLR','SSL'), data=rep('Macosko,5'), ari=rep(0,5))

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
rm(Macosko)

newX9 = X9[, sample(1:ncol(X9), 600)]
newlab = sapply(strsplit(colnames(newX9), '_'), function(x) x[1])

X9 = newX9
rm(newX9)
rm(Macosko)
rm(ind, zeros)
gc()

library(Rtsne)
#colors = rainbow(6)
#names(colors) = unique(truth9$V1)
tsne = Rtsne(t(X9), dims=2, perplexity=30, verbose=TRUE, max_iter=500)
#save(tsne, file='~/Desktop/mtmp/tsne_4800.Rdata')
#plot(tsne$Y[,1] ~ tsne$Y[,2], col=as.factor(newlab), main='true')
Y = tsne$Y

#kmeans
k        = kmeans(t(X9), numClust9, nstart=10)$cluster
k        = as.numeric(k)
plot(tsne$Y[,1] ~ tsne$Y[,2], col=as.factor(k), main='kmeans')

#pcareduce
library(pcaReduce)
pcareduce   = PCAreduce(t(X9), 10, numClust9, 'M')[[1]][,1]
plot(tsne$Y[,1] ~ tsne$Y[,2], col=as.factor(pcareduce), main='pcareduce')
#newdf$ari[2] = adj.rand.index(pcareduce, truth9$V1)

#SC3
library(gdata)
library(SC3); library(SingleCellExperiment)
sce = SingleCellExperiment(
  assays = list(
    counts = exp(as.matrix(X9))-1,
    logcounts = as.matrix(X9)
  ),
  colData = colnames(X9)
)
rowData(sce)$feature_symbol = rownames(X9)
res = sc3(sce, ks = numClust9, biology = FALSE)
sc3_export_results_xls(res, filename = paste0("Macosko.xls"))
x = read.xls("Macosko.xls")
scresult = x[,2]
plot(tsne$Y[,1] ~ tsne$Y[,2], col=as.factor(scresult), main='SC3')
#newdf$ari[3] = adj.rand.index(scresult, truth9$V1)


#SIMLR
library(SIMLR)
ss = SIMLR(X9, numClust9)
plot(tsne$Y[,1] ~ tsne$Y[,2], col=as.factor(ss$y$cluster), main='SIMLR')
ind = sort(ss$y$cluster, index.return=TRUE)$ix
heat(ss$S[ind,ind])

#newdf$ari[4] = adj.rand.index(ss$y$cluster, truth9$V1)


#SSL
set.seed(1)
R.utils::sourceDirectory("R/", modifiedOnly=FALSE)
Rcpp::sourceCpp('src/SCNoisyClustering.cpp')
P = corr_kernel(t(X9), k=30)
S = sparse_scaledlasso_c(P, 5, 0, verbose=TRUE)$S
S2 = network.diffusion(S, 30, 0.7)
#res = tsne_spectral(S2, numClust9)
res = tsne_spectral(S2, 3)
ind = sort(res, index.return=TRUE)$ix
heat(S2[ind,ind])
plot(Y[,1] ~ Y[,2], col=as.factor(res), main='SSL')
#newdf$ari[5] = adj.rand.index(res, truth9$V1)
