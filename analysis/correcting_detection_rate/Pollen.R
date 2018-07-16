library(SCNoisyClustering)
load('data/unnecessary_in_building/2_Pollen.RData')
X = as.matrix(Pollen$x)
truth = as.numeric(as.factor(Pollen$label))

heat = function(S){
  ggplot(melt(S), aes(x=X1, y=X2, fill=value)) + geom_tile() + scale_fill_gradient()
}

plot(irlba(X,2)$v[,1], irlba(X,2)$v[,2], col=truth)

logX = log(X+1)

det = colSums(X!=0) / nrow(X)
det2 = qr(det)
R = t(qr.resid(det2, t(logX)))


heat(1-cor(logX))
heat(1-cor(R))


library(Rtsne)
tsne2 = Rtsne(t(R))
tsne1 = Rtsne(t(logX))

res2 = kmeans(tsne2$Y, 11, nstart=500)$cluster
res1 = kmeans(tsne1$Y, 11, nstart=100)$cluster

adj.rand.index(res2, truth)
adj.rand.index(res1, truth)

h = hclust(as.dist(1-cor(R)))
