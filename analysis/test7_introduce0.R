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
library(CIDR)
library(pcaReduce)
library(Rtsne)


set.seed(1)
R.utils::sourceDirectory("R/", modifiedOnly=FALSE)
Rcpp::sourceCpp('src/SCNoisyClustering.cpp')
source("../SIMLR/R/tsne.R")


load('../data/3_Usoskin.RData')
X3        = Test_4_Usoskin$in_X
truth3    = Test_4_Usoskin$true_labs
numClust3 = Test_4_Usoskin$n_clust
rm(Test_4_Usoskin)

load('../data/4_Buettener.RData')
X4        = Test_1_mECS$in_X
truth4    = Test_1_mECS$true_labs; truth4 = data.frame(V1=truth4)
numClust4 = Test_1_mECS$n_clust
rm(Test_1_mECS)

load('../data/5_Yan.rda')
X5        = as.matrix(yan)
X5        = log(X5 + 1)
truth5    = as.numeric(ann$cell_type1); truth5 = data.frame(V1=truth5)
numClust5 = length(unique(truth5$V1))
rm(ann, yan)

load('../data/6_Treutlein.rda')
X6        = as.matrix(treutlein)
truth6    = as.numeric(colnames(treutlein))
ind       = sort(truth6, index.return=TRUE)$ix
X6        = X6[,ind]
X6        = log(X6 + 1)
truth6    = truth6[ind]; truth6 = data.frame(V1=truth6)
numClust6 = length(unique(truth6$V1))
rm(treutlein)

#order the data with the number of zeros for each gene
#remove the rows that have too many zeros
zerocounts = apply(X6, 1, function(x) sum(x==0))
ind = sort(zerocounts, index.return=TRUE)
X = X6[ind$ix, ]
zerocounts = ind$x
percentage_list = c(1, 0.9, 0.8, 0.7, 0.6, 0.5)
ri6 = c()
for (perc in percentage_list){
  newX = X[zerocounts < ncol(X) * perc, ]
  print(dim(newX))
  P    = mystery_kernel(t(newX), allk_input = c(20,30), sigma=c(2,1))
  S    = sparse_scaledlasso_c(P, 5, verbose=FALSE)$S
  res  = tsne_spectral(S, numClust6) #can fix numEigen argument
  ri6   = c(ri6, adj.rand.index(truth6$V1, res))
  print(ri6)
}

df = data.frame(Usoskin = ri3, Buettener=ri4, Yan=ri5, Treutlein=ri6)
df$zero_perc = percentage_list
ggplot(melt(df, id='zero_perc'), aes(x=zero_perc, y=value, col=variable)) + geom_line()


#write.table(df, 'removing_zeros.txt', col.names=TRUE, row.names=FALSE)