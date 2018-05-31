library(data.table)
#x = fread('~/Desktop/mtmp/macosko_GSE63472_P14Retina_logDGE.txt', nrows=3) #removed
dim(x)
cellnames = colnames(x)[-1]
cellnames1 = sapply(strsplit(cellnames, '_'), function(x) x[1])
types = paste0('r', 1:6)
subset = c()
for (i in 1:6){
  ind = sample(which(cellnames1==types[i]), 800)
  subset = c(subset, ind)
}

x = fread('~/Desktop/mtmp/GSE63472_P14Retina_logDGE.txt', select = c(1,subset), nrows=-1)
x = setDF(x)
label = rep(1:6, each=800)
Macosko = list(X = x, label = label)
save(Macosko, file='../data/9_Macosko.Rdata')
