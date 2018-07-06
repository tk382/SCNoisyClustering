### FINAL RESULT OF OUR METHOD ###

# We use alpha = 0.7 for network diffusion
# We use kratio = 0.05 unless the number of cells are small and we have lower bound of 10
set.seed(1)
result = matrix(0, 8, 4)
colnames(result) = c("spearman","pearson","euclidean","combined")
for (i in 2:8){
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
    X        = read.table('data/7_GSE75748_sc_cell_type_ec_cleaned.txt', header=TRUE, stringsAsFactors = FALSE)
    X        = as.matrix(X)
    truth    = data.frame(V1 = as.factor(sapply(strsplit(colnames(X), '_'), function(x) x[1])))
    truth$V1 = as.integer(truth$V1)
    numClust = 7
  }else if(i==8){
    X        = read.table('data/8_cleaned_celltype_exp2.txt', header=TRUE, stringsAsFactors = FALSE)
    X        = as.matrix(X)
    truth    = data.frame(V1 = as.factor(sapply(strsplit(colnames(X), '_'), function(x) x[1])))
    truth$V1 = as.integer(truth$V1)
    numClust = length(unique(truth$V1))
  }

  filter=TRUE
  ####################################
  res = ssl_wrapper(X, numClust, filter=filter, kerneltype="spearman",verbose=FALSE)
  result[i,1] = adj.rand.index(res$result, truth$V1)
  gc()
  res = ssl_wrapper(X, numClust, filter=filter, kerneltype="pearson",verbose=FALSE)
  result[i,2] = adj.rand.index(res$result, truth$V1)
  gc()
  res = ssl_wrapper(X, numClust, filter=filter, kerneltype="euclidean",verbose=FALSE)
  result[i,3] = adj.rand.index(res$result, truth$V1)
  gc()
  res  = ssl_wrapper(X, numClust, filter=filter, kerneltype="combined",verbose=TRUE)
  result[i,4] = adj.rand.index(res$result, truth$V1)
  gc()
  print(paste('data set',i,'is done!'))
}











load('results_references/kmeans_8_dataset.Rdata') #k
load('results_references/pcareduce_result_8_dataset.Rdata') #pcareduce
load('results_references/SC3result_clustering.Rdata') #sc3_result
load('results_references/SIMLR_result_8_dataset.Rdata') #s (these are matrices), s[[1]]$y$cluster

datanames = c('Kolod','Pollen','Usoskin','Buettener','Yan','Treutlein','Chu1','Chu2')
df = data.frame(method = rep(c('kmeans','pcareduce','SC3','SIMLR','SSL'), 8),
                data = rep(datanames, each=5),
                ari = rep(0, 40))
df[df$method=='SSL',]$ari = result
for (i in 1:8){
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
    X        = read.table('data/7_GSE75748_sc_cell_type_ec_cleaned.txt', header=TRUE, stringsAsFactors = FALSE)
    X        = as.matrix(X)
    truth    = data.frame(V1 = as.factor(sapply(strsplit(colnames(X), '_'), function(x) x[1])))
    truth$V1 = as.integer(truth$V1)
    numClust = 7
  }else if(i==8){
    X        = read.table('data/8_cleaned_celltype_exp2.txt', header=TRUE, stringsAsFactors = FALSE)
    X        = as.matrix(X)
    truth    = data.frame(V1 = as.factor(sapply(strsplit(colnames(X), '_'), function(x) x[1])))
    truth$V1 = as.integer(truth$V1)
    numClust = length(unique(truth$V1))
  }
  df[df$data==datanames[i] & df$method=='kmeans',]$ari = adj.rand.index(truth$V1, k[[i]])
  df[df$data==datanames[i] & df$method=='SC3',]$ari = adj.rand.index(truth$V1, sc3_result[[i]])
  df[df$data==datanames[i] & df$method=='SIMLR',]$ari = adj.rand.index(truth$V1, s[[i]]$y$cluster)
  df[df$data==datanames[i] & df$method=='pcareduce',]$ari = adj.rand.index(truth$V1,pcareduce[[i]])
}
ggplot(df, aes(x=data,y=ari,fill=method))+
  geom_bar(stat="identity",position="dodge")+
  coord_cartesian(ylim=c(0,1))


ggplot(newsumm, aes(x=data, y=ARI, fill=method))+
  geom_bar(stat="identity", position="dodge") +
  coord_cartesian(ylim=c(0,1))



