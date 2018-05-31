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

#### set up data ####
load('data/1_Kolod.RData')
X1        = Test_2_Kolod$in_X;
truth1    = Test_2_Kolod$true_labs
numClust1 = Test_2_Kolod$n_clust;
rm(Test_2_Kolod)

load('../data/2_Pollen.RData')
X2        = Test_3_Pollen$in_X
truth2    = Test_3_Pollen$true_labs
numClust2 = Test_3_Pollen$n_clust
rm(Test_3_Pollen)

load('data/3_Usoskin.RData')
X3        = Test_4_Usoskin$in_X
truth3    = Test_4_Usoskin$true_labs
numClust3 = Test_4_Usoskin$n_clust
rm(Test_4_Usoskin)

load('data/4_Buettener.RData')
X4        = Test_1_mECS$in_X
truth4    = Test_1_mECS$true_labs; truth4 = data.frame(V1=truth4)
numClust4 = Test_1_mECS$n_clust
rm(Test_1_mECS)

load('../data/5_Yan.rda')
X5        = as.matrix(yan)
truth5    = as.numeric(ann$cell_type1); truth5 = data.frame(V1=truth5)
numClust5 = length(unique(truth5$V1))
rm(ann, yan)

load('../data/6_Treutlein.rda')
X6        = as.matrix(treutlein)
truth6    = as.numeric(colnames(treutlein))
ind       = sort(truth6, index.return=TRUE)$ix
X6        = X6[,ind]
truth6    = truth6[ind]; truth6 = data.frame(V1=truth6)
numClust6 = length(unique(truth6$V1))
rm(treutlein)

X7        = read.table('../data/7_GSE75748_sc_cell_type_ec_cleaned.txt', header=TRUE, stringsAsFactors = FALSE)
X7        = as.matrix(X7)
truth7    = data.frame(V1 = as.factor(sapply(strsplit(colnames(X7), '_'), function(x) x[1])))
truth7$V1 = as.integer(truth7$V1)
numClust7 = 7

X8        = read.table('../data/8_cleaned_celltype_exp2.txt', header=TRUE, stringsAsFactors = FALSE)
X8        = as.matrix(X8)
truth8    = data.frame(V1 = as.factor(sapply(strsplit(colnames(X8), '_'), function(x) x[1])))
truth8$V1 = as.integer(truth8$V1)
numClust8 = length(unique(truth8$V1))

#methods to test
# kmeans
# tSNE + kmeans
# SC3
# pcaReduce
# SIMLR
# nucs

X = list(X1=X1, X2=X2, X3=X3, X4=X4, X5=X5,X6=X6,X7=X7,X8=X8)
truth = list(truth1$V1, truth2$V1, truth3$V1, truth4$V1, truth5$V1,
             truth6$V1, truth7$V1, truth8$V1)
numClust = list(numClust1, numClust2, numClust3, numClust4, numClust5,
                numClust6, numClust7, numClust8)
datanames = c('Kolod','Pollen','Usoskin','Buettener','Yan','Treutlein','Chu1','Chu2')

rm(X1,X2,X3,X4,X5,X6,X7,X8,
   truth1,truth2,truth3,truth4,truth5,truth6,truth7,truth8,
   numClust1, numClust2, numClust3, numClust4,
   numClust5, numClust6, numClust7, numClust8)

summ = data.frame(data       = rep(datanames, each=5),
                  method     = rep(c('pcaReduce','SIMLR',
                                     'sparseSL+fl','kmeans','sparseSL'), 8),
                  ARI      = rep(0, 40))


load('results_references/ours_8_dataset_clusterings.Rdata')
load('results_references/ours_8_dataset.Rdata')
load('results_references/pcareduce_result_8_dataset.Rdata')
load('results_references/SIMLR_result_8_dataset.Rdata')
load('results_references/ssl_8_dataset.Rdata')
load('results_references/ssl_8_dataset_clustering.Rdata')
load('results_references/kmeans_8_dataset.Rdata')
load('results_references/ssl_fusedlasso_nclust_plus_five_clustering.Rdata')
load('results_references/ssl_fusedlasso_nclust_plus_five.Rdata')
#SS= ssl_fl_nclust_res = list(1,2,3,4,5,6,7,8)
for (i in 1:8){
  print(i)
  tmpX         = X[[i]]
  #h            = hclust(dist(t(tmpX)))
  ind          = sort(truth[[i]], index.return=TRUE)
  # tmpX        = tmpX[, ind$ix]
  # if(i %in% c(5,6)){
  #   tmpX = log(tmpX+1)
  # }
  # gc()
  # nClust      = numClust[[i]]
  # S           = ssl[[i]]$S[ind$ix, ind$ix]
  # SS[[i]]     = fusedlasso(S, nClust+5)
  # ssl_fl_nclust_res[[i]] = tsne_spectral(SS[[i]], nClust) #can fix numEigen argument
  summ[summ$data==datanames[i] & summ$method=='kmeans', 'ARI'] =
    adj.rand.index(truth[[i]], k[[i]])
  summ[summ$data==datanames[i] & summ$method=='SIMLR', 'ARI'] =
    adj.rand.index(truth[[i]], s[[i]]$y$cluster)
  summ[summ$data==datanames[i] & summ$method=='pcaReduce', 'ARI'] =
    adj.rand.index(truth[[i]], pcareduce[[i]])
  summ[summ$data==datanames[i] & summ$method=='sparseSL', 'ARI'] =
    adj.rand.index(truth[[i]], res[[i]])
  summ[summ$data==datanames[i] & summ$method=='sparseSL+fl', 'ARI'] =
    adj.rand.index(ind$x, ssl_fl_nclust_res[[i]])
  gc()
}
ggplot(summ, aes(x=data, y=ARI, fill=method))+
  geom_bar(stat="identity", position="dodge") +
  coord_cartesian(ylim=c(0,1))

#increase the number of changepoints
