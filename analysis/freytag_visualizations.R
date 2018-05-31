#### Problem set up ####
par(mar=c(5,5,5,3))
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
set.seed(1)
R.utils::sourceDirectory("R/", modifiedOnly=FALSE)
Rcpp::sourceCpp('src/SCNoisyClustering.cpp')




#10X Gemcode - pbmc 3k
library(Matrix)
tmp = readMM('data/pbmc3k/matrix.mtx')
dim(tmp)

#pbmc 11type reference data
ref = readRDS('data/pbmc3k/all_pure_select_11types.rds')
refX = t(ref$pure_avg)
colnames(refX) = ref$pure_id
refX = scale(refX)

#take appropriate subset
Y = tmp[ref$pure_use_genes, ]
rm(tmp)


tmppca = irlba(as.matrix(Y), 2)
plot(tmppca$v[,1] ~ apply(Y, 2, function(x) sum(x==0)),
     main='relationship between coverage and pc1')
firstpc = tmppca$v[,1]

for (i in 1:ncol(Y)){
  Y[i,] = resid(lm(as.numeric(Y[i,]) ~ firstpc))
}

Y = log(Y+1)


projection = t(scale(Y)) %*% refX / nrow(Y)
projection = projection / rowSums(projection)
colnames(projection) = ref$pure_id
ggplot(melt(projection), aes(x=X2, y=X1, fill=value)) + geom_tile() + scale_fill_gradient()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


#create label on the 11 type reference set
label = as.numeric(as.factor(apply(projection, 1, function(x) colnames(projection)[which.max(x)])))

#redefine the projection data based on the label (sort them)
sorted = sort(label, index.return = TRUE)
newproj = projection[sorted$ix, ]
ggplot(melt(newproj), aes(x=X2, y=X1, fill=value)) + geom_tile() + scale_fill_gradient()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


#Y = Y[,sorted$ix]

library(Rtsne)

tsne = Rtsne(t(Y), dims=3, perplexity=30, verbose=TRUE, max_iter=500)
scatterplot3d(tsne$Y[,1], tsne$Y[,2], tsne$Y[,3], color = rainbow(11)[((label))])
png('analysis/writeup/tsne_orig_pc1out.png', width=600,height=600)
plot(tsne$Y[,1], tsne$Y[,2], col = rainbow(11)[label], main = 't-SNE of original without pc1',
     xlab='', ylab='')
dev.off()
plot_ly(x=tsne$Y[,1], y=tsne$Y[,2], z=tsne$Y[,3], color=rainbow(11)[label])

tsneproj = Rtsne(newproj, dims=3, perplexity=30, verbose=TRUE, max_iter=500)
scatterplot3d(tsneproj$Y[,1], tsneproj$Y[,2], tsneproj$Y[,3], color = rainbow(11)[as.factor(sort(label))])
png('analysis/writeup/tsne_proj_pc1out.png', width=600,height=600)
plot(tsneproj$Y[,1], tsneproj$Y[,2], col= rainbow(11)[sort(label)], main='t-SNE of projection onto reference without pc1',
     xlab='', ylab='')
dev.off()w
label2 = sort(label)
plot_ly(x=tsneproj$Y[,1], y=tsneproj$Y[,2], z=tsneproj$Y[,3],
       color=rainbow(11)[sort(label)])


zeros = apply(Y, 2, function(x) sum(x==0))
sv = irlba(Y, 5)

png('analysis/writeup/pc1_zeros.png', width=600, height=600)
plot(sv$v[,1] ~ zeros,
     xlab='# of zeros for each cell',
     ylab='pc1')
dev.off()

png('analysis/writeup/pca_orig_pc1out.png', width=600,height=600)
plot(sv$v[,2] ~ sv$v[,3], col = rainbow(11)[label], main='pca of original without pc1',
     xlab = 'pc2', ylab = 'pc1')
dev.off()

sv2 = irlba(t(newproj), 5)
png('analysis/writeup/pca_proj_pc1out.png', width=600,height=600)
plot(sv2$v[,1] ~ sv2$v[,2], col=rainbow(11)[sort(label)], main='pca of projection onto 11 reference',
     xlab = 'pc2', ylab = 'pc1')
dev.off()





######projection onto the 109 cells ##########
load('data/sysdata.rda')
Xgenenames = read.table('~/Desktop/mtmp/pbmc_3k_filtered_gene_bc_matrices/hg19/genes.tsv')
intersection = (which(rownames(sysdata$GlobalPanel[[1]]) %in% Xgenenames$V2))
panel = sysdata$GlobalPanel[[1]][intersection, ]
ind = match(rownames(panel), Xgenenames$V2)
X = readMM('~/Desktop/mtmp/pbmc_3k_filtered_gene_bc_matrices/hg19/matrix.mtx')
X = X[ind, ]

projection = t(scale(X)) %*% scale(as.matrix(panel)) / nrow(X)
projection = projection / rowSums(projection)
png('analysis/writeup/heatmap_projection_onto_cell_atlas.png', width=1000,height=500)
ggplot(melt(projection), aes(x=X2, y=X1, fill=value)) + geom_tile() + scale_fill_gradient(high = 'white', low = 'navy')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ylab('')+xlab('')+
  ggtitle("Projection onto Primary Cell Atlas")
dev.off()



bigt = Rtsne(projection, dims=5, perplexity=30, verbose=TRUE, max_iter=500)
png('analysis/writeup/tsne_proj_cell_atlas.png', width=600, height=600)
plot(bigt$Y[,1], bigt$Y[,2], col = rainbow(11)[label],
     main = 't-SNE of Projection onto Primary Cell Atlas', ylab='', xlab='')
#scatterplot3d(bigt$Y[,4], bigt$Y[,3], bigt$Y[,1], color=rainbow(11)[label])
dev.off()



intersection = (which(rownames(sysdata$GlobalPanel[[2]]) %in% Xgenenames$V2))
panel = sysdata$GlobalPanel[[2]][intersection, ]
ind = match(rownames(panel), Xgenenames$V2)
X = readMM('~/Desktop/mtmp/pbmc_3k_filtered_gene_bc_matrices/hg19/matrix.mtx')
X = X[ind, ]; gc()

projection = t(scale(X)) %*% scale(as.matrix(panel)) / nrow(X)
projection = projection / rowSums(projection)
png('analysis/writeup/heatmap_projection_onto_gene_atlas.png', width=1000,height=500)
ggplot(melt(projection), aes(x=X2, y=X1, fill=value)) + geom_tile() + scale_fill_gradient(high = 'white', low = 'navy')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+xlab('')+ylab('cells')+
  ggtitle('Projection onto GNF1H Gene Atlas')
dev.off()
bigt2 = Rtsne(projection, dims=5, perplexity=30, verbose=TRUE, max_iter=500)
png('analysis/writeup/tsne_proj_gene_atlas.png', width=600, height=600)
plot(bigt2$Y[,1], bigt2$Y[,2], col = rainbow(11)[label],
     main = 't-SNE of Projection onto GNF1H Gene Atlas', xlab='', ylab='')
dev.off()



scatterplot3d(bigt$Y[,3], bigt$Y[,4], bigt$Y[,5], color=rainbow(11)[label])
