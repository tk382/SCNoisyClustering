### FINAL RESULT OF OUR METHOD ###

# We use alpha = 0.7 for network diffusion
# We use kratio = 0.05 unless the number of cells are small and we have lower bound of 10
result = c()
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
  ssl_wrapper(X, numClust)
  # zeros = apply(X, 1, function(x) sum(x==0))
  # remove = which(zeros > ncol(X)*0.95)
  # if(length(remove)>0){X = X[-remove,]}
  # k     = max(round(ncol(X) * 0.05), 10)
  # allk  = seq(10,30,by=5); allsigma = seq(2,1,by=-0.1)
  # P = corr_kernel(t(X), k, allk, allsigma)
  # res = sparse_scaledlasso_c(P, 5, 1, verbose=TRUE)
  # S = network.diffusion(res$S, k, alpha=0.7)
  # estimates = tsne_spectral(S, numClust, 5)
  ari = adj.rand.index(estimates, truth$V1)
  result = c(result, ari)
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



