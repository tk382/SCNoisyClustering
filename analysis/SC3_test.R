library(SC3)
library(SingleCellExperiment)
library(scater)
library(gdata)

set.seed(1)
R.utils::sourceDirectory("R/", modifiedOnly=FALSE)
Rcpp::sourceCpp('src/SCNoisyClustering.cpp')
source("../benchmark/SIMLR/R/tsne.R")

#### set up data ####
load('data/1_Kolod.RData')
X1        = Test_2_Kolod$in_X;
truth1    = Test_2_Kolod$true_labs
numClust1 = Test_2_Kolod$n_clust;
rm(Test_2_Kolod)

load('data/2_Pollen.RData')
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

load('data/5_Yan.rda')
X5        = as.matrix(yan)
truth5    = as.numeric(ann$cell_type1); truth5 = data.frame(V1=truth5)
numClust5 = length(unique(truth5$V1))
rm(ann, yan)

load('data/6_Treutlein.rda')
X6        = as.matrix(treutlein)
truth6    = as.numeric(colnames(treutlein))
ind       = sort(truth6, index.return=TRUE)$ix
X6        = X6[,ind]
truth6    = truth6[ind]; truth6 = data.frame(V1=truth6)
numClust6 = length(unique(truth6$V1))
rm(treutlein)

X7        = read.table('data/7_GSE75748_sc_cell_type_ec_cleaned.txt', header=TRUE, stringsAsFactors = FALSE)
X7        = as.matrix(X7)
truth7    = data.frame(V1 = as.factor(sapply(strsplit(colnames(X7), '_'), function(x) x[1])))
truth7$V1 = as.integer(truth7$V1)
numClust7 = 7

X8        = read.table('data/8_cleaned_celltype_exp2.txt', header=TRUE, stringsAsFactors = FALSE)
X8        = as.matrix(X8)
truth8    = data.frame(V1 = as.factor(sapply(strsplit(colnames(X8), '_'), function(x) x[1])))
truth8$V1 = as.integer(truth8$V1)
numClust8 = length(unique(truth8$V1))

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


res = result = list(1,2,3,4,5,6,7,8)
for (i in 7:8){
  if(i %in% 5:6){
    colnames(X[[i]]) = paste0('C',1:ncol(X[[i]]))
    rownames(X[[i]]) = paste0('R',1:nrow(X[[i]]))
    sce = SingleCellExperiment(
      assays = list(
        counts = as.matrix(X[[i]]),
        logcounts = log2(as.matrix(X[[i]]) + 1)
      ),
      colData = colnames(X[[i]])
    )
  }else{
    colnames(X[[i]]) = paste0('C',1:ncol(X[[i]]))
    rownames(X[[i]]) = paste0('R',1:nrow(X[[i]]))
    sce = SingleCellExperiment(
      assays = list(
        counts = exp(as.matrix(X[[i]]))-1,
        logcounts = as.matrix(X[[i]])
      ),
      colData = colnames(X[[i]])
    )
  }

  rowData(sce)$feature_symbol = rownames(X[[i]])
  res[[i]] = sc3(sce, ks = numClust[[i]], biology = FALSE)
  sc3_export_results_xls(res[[i]], filename = paste0(datanames[i], ".xls"))
  x = read.xls(paste0(datanames[i], ".xls"))
  result[[i]] = x[,2]
  print(adj.rand.index(truth[[i]], result[[i]]))
}
save(sc3_res, file = "results_references/SC3result.Rdata")
save(sc3_result, file = "results_references/SC3result_clustering.Rdata")

