getwd()
# load('data/unnecessary_in_building/2_Pollen.RData')
# X        = as.matrix(Pollen$x)
# zeros = colSums(X==0)
# ind = which(zeros < 3500)
# X = X[,ind]
# truth    = Pollen$label
# truth = truth[ind]
# numClust = 11
# rm(Pollen)
# X = X[,-c(67,68)]
# truth = truth[-c(67,68)]
# truth = as.numeric(as.factor(truth))


load('data/unnecessary_in_building/Zeisel.RData')
X = Zeisel$X
truth = as.numeric(as.factor(Zeisel$label))
logX = log(X+1)

det = colSums(X!=0) / nrow(X)
# det2 = qr(det)
# R = t(qr.resid(det2, t(logX)))

pc1 = irlba(logX,1)$v[,1]
pc11 = qr(pc1)
R = t(qr.resid(pc11, t(logX)))

library(irlba)
plot(irlba(logX,1)$v[,1], det)
plot(irlba(R,1)$v[,1], det)



#make a subset of the cells
label = c()
for (i in 1:9){
  if(i==8){
    label = c(label, which(truth==i))
  }else{
    ind = which(truth==i)
    label = c(label, sample(ind, 30))
  }
}

logX = logX[, label]
truth = truth[label]

R = R[, label]

#clustering on R
zeros = rowSums(logX==0) / ncol(logX)
ind = which(zeros >= 0.9)
R = R[-ind,]

logX = genefilter(logX)

#build kernel
diff1 = 1-cor(as.matrix(R), method = "pearson")
diff3 = 1-cor(as.matrix(R), method = "spearman")
P1 = corr_kernel_c(t(R), diff1, c(15,20,25), seq(1,2,by=0.2), 20)
P2 = dist_kernel_c(t(R), c(15,20,25), seq(1,2,by=0.2), 20)
P3 = rank_kernel_c(t(R), diff3, c(15,20,25), seq(1,2,by=0.2), 20)
P  = abind(P1, P2, P3)
rm(P1, P2, P3)



diff1 = 1-cor(as.matrix(logX), method = "pearson")
diff3 = 1-cor(as.matrix(logX), method = "spearman")
P1 = corr_kernel_c(t(logX), diff1, c(15,20,25), seq(1,2,by=0.2), 20)
P2 = dist_kernel_c(t(logX), c(15,20,25), seq(1,2,by=0.2), 20)
P3 = rank_kernel_c(t(logX), diff3, c(15,20,25), seq(1,2,by=0.2), 20)
PP  = abind(P1, P2, P3)
rm(P1, P2, P3)



Y1 = Rtsne(P[,,19])$Y
Y2 = Rtsne(PP[,,19])$Y


V1 = irlba(P[,,19],2)$v
V2 = irlba(PP[,,19],2)$v

plot(V1, col=rainbow(11)[truth])
plot(V2[,2], -V2[,1], col=rainbow(11)[truth])

