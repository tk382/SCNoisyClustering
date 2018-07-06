result = matrix(0,8,4)
colnames(result) = c('spearman','pearson','euclid','combined')
rownames(result) = c('kolod','pollen','usoskin','buettner','yan','treutlein','chu1','chu2')

result2 = result


for (i in  2:8){
  set.seed(1)
  if(i==2){
    load('data/2_Pollen.RData')
    X        = Pollen$x
    zeros = colSums(X==0)
    ind = which(zeros>=5300)
    X = X[,-ind]
    truth    = Pollen$label
    truth = truth[-ind]
    numClust = 11
    rm(Pollen)
  }
  else if(i==3){
    load('data/3_Usoskin.RData')
    X        = Usoskin$X
    #X        = log(X+1)
    truth    = as.character(Usoskin$lab1)
    numClust = 4
    rm(Usoskin)
  }
  else if(i==4){
    load('data/4_Buettner.RData')
    X        = Buettner$X
    #X        = log(X+1)
    truth    = as.character(Buettner$label)
    numClust = 3
    rm(Buettner)
    ind1 = which(colSums(X==0) > 35000)
    X = X[,-ind1]; truth = truth[-ind1]
    v2 = irlba(X,2)$v[,2]
    ind2 = which(abs(v2-mean(v2)) > iqr(v2) * 1.5)
    X = X[,-ind2]; truth = truth[-ind2]
  }
  else if(i==5){
    load('data/5_Yan.rda')
    X        = as.matrix(yan)
    #X        = log(X+1)
    truth    = as.character(ann$cell_type1)
    numClust = 6
    rm(ann, yan)
  }
  else if(i==6){
    load('data/6_Treutlein.rda')
    X        = as.matrix(treutlein)
    truth    = as.numeric(colnames(treutlein))
    ind      = sort(truth, index.return=TRUE)$ix
    X        = X[,ind]
    #X        = log(X + 1)
    truth    = truth[ind]
    numClust = length(unique(truth))
    rm(treutlein)
  }
  else if(i==7){
    load('data/7_Chu_celltype.Rdata')
    X        = Chu_celltype$X
    truth    = Chu_celltype$label
    numClust = 7
  }
  else{
    load('data/8_Chu_timecourse.Rdata')
    X        = Chu_timecourse$X
    truth    = Chu_timecourse$label
    numClust = length(unique(truth))
  }


  filter=TRUE

  ####################################
  heat = function(S){
    ggplot(melt(S), aes(x=X1, y=X2, fill=value)) + geom_tile() + scale_fill_gradient()
  }
  set.seed(1)
  truth = as.numeric(as.factor(truth))
  # res = SLSL(X=X, numClust = numClust, kernel_type="spearman",verbose=T)
  # S_spearman = res$S
  # result[i,1] = adj.rand.index(res$result, truth)
  # gc()
  # res = SLSL(X, numClust =numClust, kernel_type="pearson",verbose=F)
  # S_pearson = res$S
  # result[i,2] = adj.rand.index(res$result, truth)
  # gc()
  # res = SLSL(X, numClust =numClust, kernel_type="euclidean",verbose=F)
  # S_euclid = res$S
  # result[i,3] = adj.rand.index(res$result, truth)
  # gc()
  # res = SLSL(X, numClust =numClust, kernel_type="combined",verbose=T)
  # S_combined = res$S
  # result[i,4] = adj.rand.index(res$result, truth)
  # gc()

#
#   truth = as.numeric(as.factor(truth))
#   res = SLSL(X=X, numClust = numClust, kernel_type="spearman",verbose=F)
#   result2[i,1] = adj.rand.index(res$result, truth)
#   gc()
#   res = SLSL(X, numClust =numClust, kernel_type="pearson",verbose=F)
#   result2[i,2] = adj.rand.index(res$result, truth)
#   gc()
#   res = SLSL(X, numClust =numClust, kernel_type="euclidean",verbose=F)
#   result2[i,3] = adj.rand.index(res$result, truth)
#   gc()
  res = SLSL(X, numClust =numClust, gamma=1, kernel_type="combined",verbose=T)
  result2[i,4] = adj.rand.index(res$result, truth)
  gc()
}

#
# plot_ly(x=res_tsne$vis[,1], y=res_tsne$vis[,2], z=res_tsne$vis[,3], color=truth)
# plot_ly(x=res_pca$vis[,1],  y=res_pca$vis[,2],  z=res_pca$vis[,3], color=truth,
#         main = 'Usoskin using PCA')
