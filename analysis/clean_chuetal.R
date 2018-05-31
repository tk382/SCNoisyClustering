library(matrixStats)
library(irlba)
library(SIMLR)
library(ggplot2)
library(KRLS)
library(dplyr)
library(reshape)
library(caret)
library(fossil)
library(gridExtra)
library(pracma)
library(igraph)
library(parallel)
library(scatterplot3d)

ct = read.csv('../data/GSE75748_sc_cell_type_ec.csv')
cells = colnames(ct)
spl = strsplit(cells, '_')
spl2 = sapply(spl, function(x) x[2])
spl2 = strsplit(spl2, '[.]')
celltypes = sapply(spl, function(x) x[1])
exp = sapply(spl2, function(x) x[1])
celltypes = celltypes[-1]
exp = exp[-1]
length(unique(celltypes))

#exp1
result = data.frame(X = ct[,1])
picked_celltypes = "X"
for (i in 1:length(unique(celltypes))){
  ind = which(celltypes == unique(celltypes)[i])
  exps = unique(exp[ind])
  finalwhere = ind[which(exp[ind]==unique(exp[ind])[1])]
  result = cbind(result, ct[, (finalwhere+1)])
  picked_celltypes = c(picked_celltypes, celltypes[finalwhere])
}

#exp_j
j=3
resultj = data.frame(X = ct[,1])
picked_celltypes = "X"
for (i in 1:length(unique(celltypes))){
  ind = which(celltypes == unique(celltypes)[i])
  exps = unique(exp[ind])
  finalwhere = ind[which(exp[ind]==unique(exp[ind])[j])]
  resultj = cbind(resultj, ct[, (finalwhere+1)])
  picked_celltypes = c(picked_celltypes, celltypes[finalwhere])
}


tempX = resultj[,-1]
logX = log(tempX+1)

vars = apply(logX, 1, var)
ind = which(vars > 1)
X = logX[ind, ]

write.table(X, file='cleaned_celltype_exp3.txt', col.names=TRUE, row.names=FALSE, quote = FALSE)

