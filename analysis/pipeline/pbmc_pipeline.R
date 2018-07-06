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

library(diceR)


set.seed(1)
R.utils::sourceDirectory("R/", modifiedOnly=FALSE)
Rcpp::sourceCpp('src/SCNoisyClustering.cpp')




#10X Gemcode - pbmc 3k
library(Matrix)
orig = readMM('data/pbmc3k/matrix.mtx')
genenames = read.table('data/pbmc3k/genes.tsv')
dim(orig)

#take log
lorig = log(orig + 1)

#first pc
lorig_pc1 = irlba(lorig, 1)$v[,1]

#zero counts
zeros = apply(lorig, 2, function(x) sum(x==0))

#check dependency
plot(lorig_pc1 ~ zeros, xlab='zero counts', ylab='pc1')

#regress out and get residuals
x = zeros - mean(zeros)
lorig_rs = t(t(scale(lorig)) - x %*% (t(x) %*% t(scale(lorig))/(sum(x^2))))

lrs_pc1 = irlba(lorig_rs, 1)$v[,1]
lrs_zeros = apply(lorig, 2, function(x) sum(x==0))
plot(lrs_pc1 ~ lrs_zeros, xlab='detection rate',
     ylab = 'pc1',
     main='log count matrix after removing pc1')

#split data
splits = ceil(ncol(lorig)/300)
indlist = split(1:2700, sample(rep(1:9, 300)))
X = list()
for (i in 1:9){
  X[[i]] = lorig_rs[,indlist[[i]]]
}




load('data/sysdata.rda')
primary_cell_atlas = sysdata$GlobalPanel[[1]]
gene_atlas = sysdata$GlobalPanel[[2]]
rm(sysdata)

groupslist = combn(1:9, 2)
N = ncol(groupslist)
cluster_result_geneatlas = matrix(NA, 2700, N)
cluster_result_cellatlas = matrix(NA, 2700, N)

for (i in 2:N){
  groups = groupslist[,i]
  x = cbind(X[[groups[1]]], X[[groups[2]]])
  rownames(x) = genenames$V2
  int1 = base::intersect(rownames(primary_cell_atlas), rownames(x)) #3680
  int2 = base::intersect(rownames(gene_atlas), rownames(x)) #4004

  ind1 = match(int1, rownames(x))
  x1 = x[ind1, ]
  ind11 = match(int1, rownames(primary_cell_atlas))
  primary_cell_atlas = primary_cell_atlas[ind11, ]

  ind2 = match(int2, rownames(x))
  x2 = x[ind2, ]
  ind22 = match(int2, rownames(gene_atlas))
  gene_atlas = gene_atlas[ind22, ]

  projection1 = t(t(x1) %*% as.matrix(primary_cell_atlas))/nrow(x1)
  colnames(projection1) = factor(colnames(projection1), levels = colnames(projection1))
  projection1 = scale(projection1^4)
  #heat(scale(projection1^4),'')
  res = ssl_wrapper(projection1)
  cluster_result_geneatlas[c(indlist[[groups[1]]], indlist[[groups[2]]]),i] = res$result


  projection2 = t(t(x2) %*% as.matrix(gene_atlas))/nrow(x2)
  colnames(projection2) = factor(colnames(projection2), levels = colnames(projection2))
  projection2 = scale(projection2^4)
  res = ssl_wrapper(projection2)
  cluster_result_cellatlas[c(indlist[[groups[1]]], indlist[[groups[2]]]),i] = res$result

  gc()
}

final = array(0, dim=c(nrow(cluster_result_cellatlas), ncol(cluster_result_cellatlas),1,1))
final[,,1,1]= cluster_result_cellatlas
k=8
dimnames(final)[[1]] = paste0('R',1:nrow(cluster_result_cellatlas))
dimnames(final)[[2]] = paste0('C',1:ncol(cluster_result_cellatlas))
dimnames(final)[[3]] = paste0("k")
dimnames(final)[[4]] = as.character(k)

final_cell = CSPA(final, 8)

final = array(0, dim=c(nrow(cluster_result_geneatlas), ncol(cluster_result_geneatlas),1,1))
final[,,1,1]= cluster_result_geneatlas
k=4
dimnames(final)[[1]] = paste0('R',1:nrow(cluster_result_geneatlas))
dimnames(final)[[2]] = paste0('C',1:ncol(cluster_result_geneatlas))
dimnames(final)[[3]] = paste0("k")
dimnames(final)[[4]] = as.character(k)

final_gene = CSPA(final, 4)

dim(projection1)
dim(projection2)

library(Rtsne)

cell_tsne = Rtsne(t(projection1), dims=3, perplexity=20, verbose=TRUE, max_iter=1500)
gene_tsne = Rtsne(t(projection2), dims=3, perplexity=20, verbose=TRUE, max_iter=1500)

library(plotly)
plot_ly(x=cell_tsne$Y[,1], y=cell_tsne$Y[,2], z=cell_tsne$Y[,3], color = final_cell)
plot_ly(x=gene_tsne$Y[,1], y=gene_tsne$Y[,2], z=gene_tsne$Y[,3], color = final_gene)
