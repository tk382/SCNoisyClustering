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
library(pcaReduce)
library(Rtsne)


set.seed(1)
R.utils::sourceDirectory("R/", modifiedOnly=FALSE)
Rcpp::sourceCpp('src/SCNoisyClustering.cpp')


i = 2:6
alpha = seq(0.5,0.9,by=0.1)
dist = c("dist","corr")
kratio = seq(0.05, 0.2, by=0.05)
df = expand.grid(i = i, alpha=alpha, dist=dist, kratio=kratio)
df$ri = rep(0, nrow(df))

for (rowcount in 1:nrow(df)){
  i = df$i[rowcount]
  dist = df$dist[rowcount]
  kratio = df$kratio[rowcount]
  alpha = df$alpha[rowcount]
  if (i==1){
    load('data/1_Kolod.RData')
    X        = Test_2_Kolod$in_X;
    truth    = Test_2_Kolod$true_labs
    numClust = Test_2_Kolod$n_clust;
    rm(Test_2_Kolod)
  }else if(i==2){
    load('data/2_Pollen.RData')
    X        = Test_3_Pollen$in_X
    truth    = Test_3_Pollen$true_labs
    ind      = sort(truth$V1, index.return=TRUE)$ix
    X        = X[,ind]
    truth$V1 = truth$V1[ind]
    numClust = Test_3_Pollen$n_clust
    rm(Test_3_Pollen)
  }else if(i==3){
    load('data/3_Usoskin.RData')
    X        = Test_4_Usoskin$in_X
    truth    = Test_4_Usoskin$true_labs
    numClust = Test_4_Usoskin$n_clust
    rm(Test_4_Usoskin)
  }else if(i==4){
    load('data/4_Buettener.RData')
    X        = Test_1_mECS$in_X
    truth    = Test_1_mECS$true_labs; truth = data.frame(V1=truth)
    numClust = Test_1_mECS$n_clust
    rm(Test_1_mECS)
  }else if(i==5){
    load('data/5_Yan.rda')
    X        = as.matrix(yan)
    X        = log(X+1)
    truth    = as.numeric(ann$cell_type1); truth = data.frame(V1=truth)
    numClust = length(unique(truth$V1))
    rm(ann, yan)
  }else if(i==6){
    load('data/6_Treutlein.rda')
    X        = as.matrix(treutlein)
    truth    = as.numeric(colnames(treutlein))
    ind      = sort(truth, index.return=TRUE)$ix
    X        = X[,ind]
    X        = log(X + 1)
    truth    = truth[ind]; truth = data.frame(V1=truth)
    numClust = length(unique(truth$V1))
    rm(treutlein)
  }else if(i==7){
    X        = read.table('../data/7_GSE75748_sc_cell_type_ec_cleaned.txt', header=TRUE, stringsAsFactors = FALSE)
    X        = as.matrix(X)
    truth    = data.frame(V1 = as.factor(sapply(strsplit(colnames(X7), '_'), function(x) x[1])))
    truth$V1 = as.integer(truth$V1)
    numClust = 7
  }else if(i==8){
    X        = read.table('../data/8_cleaned_celltype_exp2.txt', header=TRUE, stringsAsFactors = FALSE)
    X        = as.matrix(X)
    truth    = data.frame(V1 = as.factor(sapply(strsplit(colnames(X), '_'), function(x) x[1])))
    truth$V1 = as.integer(truth$V1)
    numClust = length(unique(truth$V1))
  }

  #gene filter
  zeros = apply(X, 1, function(x) sum(x==0))
  remove = which(zeros > ncol(X)*0.95)
  if(length(remove)>0){X = X[-remove,]}
  k     = round(ncol(X) * kratio)
  allk  = seq(10,30,by=5); allsigma = seq(2,1,by=-0.1)
  if(dist == "dist"){
    P = mystery_kernel(t(X), k, allk, allsigma)
  }else if (dist == "corr"){
    P = corr_kernel(t(X), k, allk, allsigma)
  }
  res = sparse_scaledlasso_c(P, 5, 1, verbose=TRUE)
  S = network.diffusion(res$S, k, alpha)
  estimates = tsne_spectral(S, numClust, 5)
  ri = adj.rand.index(estimates, truth$V1)
  df$ri[rowcount] = ri
}

write.table(df, 'finetuning_parameters.txt', col.names=TRUE, row.names=FALSE, quote=FALSE)


df = read.table('analysis/finetuning_parameters.txt', header=TRUE, stringsAsFactors = FALSE)
library(ggplot2)
library(dplyr)
df = as.tbl(df)
df %>% dplyr::filter(i==2)
  ggplot(df %>% dplyr::filter(i==6,alpha==0.8), aes(x=kratio, y=ri, col=dist))+ geom_line()

hist(df[df$dist=='corr', ]$ri)
hist(df[df$dist=='dist', ]$ri)


dfdist = df[df$dist=='dist',]
dfcorr = df[df$dist=='corr',]

dfdist = dfdist[sort(dfdist$i,index.return=TRUE)$ix, ]
dfcorr = dfcorr[sort(dfcorr$i,index.return=TRUE)$ix, ]
plot(dfdist$ri[dfdist$alpha==0.7], type='l', col='blue', ylim=c(0,1))
lines(dfcorr$ri[dfdist$alpha==0.7], type='l', col='red', lty=2)



n = c(249,622,182,90,80)
#i=2 : kratio=0.05, alpha doesn't matter
#i=3 : alpha = 0.6, 0.7, 0.8, 0.7 is the best, 0.8 is the next
#    : kratio is 0.05
#i=4 : alpha = 0.8 is the best, 0.5, 6, 7 are very close. NOT 0.9
#i=5 : 0.5, 6, 7, 8 all ok
#i=6 : kratio = 0.15, alpha = 0.7 is the best. 0.5, 6 also okay.


#choose k = max(0.05 * ncol(X) , 10)
#choose alpha = 0.7.

dfcorr[c(3, 23, 43, 68, 93), ]
