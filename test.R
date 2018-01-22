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

set.seed(1)
R.utils::sourceDirectory("R/", modifiedOnly=FALSE)
#Rcpp::sourceCpp('src/temp.cpp')

#which data set to test?
data_index = 1; data_name = 'mECS'
# data_index = 2; data_name = 'Kolod' #large
# data_index = 3; data_name = 'Pollen'
# data_index = 4; data_name = 'Usoskin'
# data_index = 5; data_name = 'Zeisel' #large

#load data and define X, true labels, and number of clusters
datname = paste0('Test_',data_index, '_',data_name)
load(paste0('../data/',datname,'.Rdata'))
assign("dat", get(datname))
X = dat$in_X
truth = dat$true_labs
numClust = dat$n_clust
rm(data_index, data_name, datname)
rm(Test_1_mECS)

K = multiple_kernel_new(t(X), 1)
P = array(0,dim=c(ncol(X),ncol(X),22))
for (i in 1:22){
  P[,,i] = as.matrix(K[[i]])
}
n = 22

D = as.matrix(dist(t(X)))

tempK = constructKernel(t(X), 5)
diag(tempK) = colMeans(tempK)

theirs = ggplot(melt(as.matrix(K[[1]])), aes(x=X1, y=X2, fill=value)) + geom_tile() + scale_fill_gradient()
ours = ggplot(melt(tempK), aes(x=X1, y=X2, fill=value)) + geom_tile() + scale_fill_gradient()
grid.arrange(theirs, ours, nrow=1)


skip = round(nrow(X)/3)
X1 = X[1:skip, ]
X2 = X[(skip+1) : (2*skip),]
X3 = X[(2*skip+1) : nrow(X), ]
K1 = multiple_kernel(t(X1), 1)[[21]]
K2 = multiple_kernel(t(X2), 1)[[21]]
K3 = multiple_kernel(t(X3), 1)[[21]]
K = array(0, dim=c(ncol(X), ncol(X), 3))
K[,,1] = as.matrix(K1)
K[,,2] = as.matrix(K2)
K[,,3] = as.matrix(K3)
rm(K1, K2, K3)

#T is similarity matrices created from different features
T = array(0, dim = c(ncol(X), ncol(X), 3))
for (i in 1:3){
  Ki = K[,,i]
  D = diag(1/rowSums(Ki))
  Ti = D %*% Ki
  T[,,i] = Ti
}
D = diag(1/rowSums(Ki))
Ti = D %*% Ki
T[,,i] = Ti
rm(K)

#run RMSC
res_mv = RMSC_multiview(T, 0.001, 3, mu=0.001, verbose=TRUE)
S_mv = res_mv$P
E_mv = res_mv$E
res_mv = baseline_spectral(S_mv, numClust=numClust)
rand.index(truth$V1, res_mv)
compare(truth$V1, res_mv, method='nmi')
table(truth$V1, res_mv)


res_fr = RMSC_fixedrank(P, 3, mu=1e-3, rho = 1.9, max_iter=100, eps=1e-9, verbose=TRUE)
S_fr = res_fr$S
f_fr = res_fr$f
sigma_fr = res_fr$sigma
res_fr = baseline_spectral(S_fr, numClust=numClust)
rand.index(truth$V1, res_fr)
compare(truth$V1, res_fr, method='nmi')
table(truth$V1, res_fr)
ggplot(melt(S_fr), aes(x=X1, y=X2, fill=value)) + geom_tile() + scale_fill_gradient()
plot(sigma)
plot(w)

res_wt = RMSC_weighted(P, 1e-5, 1e-5, mu=1e-3, rho = 1.9, max_iter=100, eps=1e-9, verbose=TRUE)
S_wt = res_wt$S
E_wt = res_wt$E
f_wt = res_wt$f
sigma_wt = res_wt$sigma
res_wt = baseline_spectral(S_wt, numClust=numClust)
rand.index(truth$V1[random_subsample], res_wt)
compare(truth$V1[random_subsample], res_wt, method='nmi')
table(truth$V1[random_subsample], res_wt)

res_w = RMSC_with_w_without_sigma(P,numClust)
S_w = res_w$S
w_w = res_w$w
res_w2 = baseline_spectral(S_w, numClust=numClust)
rand.index(truth$V1, res_w2)
compare(truth$V1, res_w2, method='nmi')
table(truth$V1, res_w2)
ggplot(melt(S_w), aes(x=X1, y=X2, fill=value)) + geom_tile() + scale_color_gradient()

res_ws = RMSC_with_w(P, numClust)


res_sl = RMSC_without_scaled_lasso(P, numClust)
S_sl = res_sl$S
w_sl = res_sl$w
res_sl = baseline_spectral(S_sl, numClust=numClust)
rand.index(truth$V1, res_sl)
compare(truth$V1, res_sl, method='nmi')
table(truth$V1, res_sl)
ggplot(melt(S_sl), aes(x=X1, y=X2, fill=value)) + geom_tile() + scale_color_gradient()



#make clusters based on P
scatterplot3d(x = svd(S_mv)$v[,1], y=svd(S_mv)$v[,2], z=svd(S_mv)$v[,3])
plot(svd(S_fr)$v[,1:3])
#table(truth$V1, res1)

library(SIMLR)
library(Matrix)
library(parallel)
library(igraph)
library(grDevices)
source("../SIMLR/R/SIMLR.R")
source("../SIMLR/R/compute.multiple.kernel.R")
source("../SIMLR/R/network.diffusion.R")
source("../SIMLR/R/utils.simlr.R")
source("../SIMLR/R/tsne.R")
dyn.load("../SIMLR/R/projsplx_R.so")
load(file="../SIMLR/data/Test_1_mECS.RData")
res_example1 = SIMLR(X,3)
S_simlr = res_example1$S
w_simlr = res_example1$alphaK
scatterplot3d(x = svd(S_simlr)$v[,1], y=svd(S_simlr)$v[,2], z=svd(S_simlr)$v[,3])
plot(svd(D)$v[,1:3])
res_simlr = baseline_spectral(S_simlr, 3)
rand.index(truth$V1[random_subsample], res_simlr)
compare(truth$V1[random_subsample], res_simlr, method='nmi')
table(truth$V1[random_subsample], res_simlr)
#RI
#Buettner : 0.906988
#Kolod    : 0.8046683
#Pollen   : 0.9802112
#Usoskin  : 0.7171298

#NMI
#Buettner : 0.7752997
#Kolod    : 0.5542772
#Pollen   : 0.9144992
#Usoskin  : 0.3259396

E1 = E[,,1]
diag(E1) = colMeans(E1)

E2 = E[,,2]
diag(E2) = colMeans(E2)

ggplot(melt(D), aes(x=X1, y=X2, fill=value)) + geom_tile() +
  scale_color_gradient()
ggplot(melt(E[,,1]), aes(x=X1, y=X2, fill=value)) + geom_tile() +
  scale_color_gradient()
ggplot(melt(E[,,2]), aes(x=X1, y=X2, fill=value)) + geom_tile() +
  scale_color_gradient()
ggplot(melt(t(X) %*% X), aes(x=X1, y=X2, fill=value)) + geom_tile() +
  scale_color_gradient()+ggtitle('cross_prod')

pc = svd(X)
kk = kmeans(pc$v[, 1:3], 3)

E2 = E[,,2]
diag(E2) = colMeans(E2)
ggplot(melt(S), aes(x=X1, y=X2, fill=value)) +
  geom_tile() +
  scale_fill_gradient(high='indianred', low='antiquewhite1')


##sensitivity_n
ri_n = rep(0,9)
result_n = matrix(0,ncol(X), 9)
reshuffle = sample(1:nrow(X))
for (n in 2:10){
  T = make_T_matrix(X, reshuffle, n, sigma=5)
  tic(); P = RMSC_c(T, 0.05, mu=1e-3)$P; toc()
  res = baseline_spectral(P, numClust=numClust)
  result_n[,n-1] = res
  ri_n[n-1] = rand.index(truth$V1, res)
}


##sensitivity_lambda
ri_l = rep(0,9)
result_l = matrix(0,ncol(X), 9)
lambdavec = seq(0.001, 1, length=9)
n=3
reshuffle = sample(1:nrow(X))
for (i in 1:9){
  T = make_T_matrix(X, reshuffle, n, sigma=5)
  tic(); P = RMSC_c(T, lambdavec[i], mu=1e-3)$P; toc()
  res = baseline_spectral(P, numClust=numClust)
  result_l[,i] = res
  ri_l[i] = rand.index(truth$V1, res)
}



##sensitivity_kernel
ri_k = rep(0,9)
result_k = matrix(0,ncol(X), 9)
sigmavec = 1:9
n=3
reshuffle = sample(1:nrow(X))
for (i in 1:9){
  sigma = sigmavec[i]
  T = make_T_matrix(X, reshuffle, n, sigma)
  tic(); P = RMSC_c(T, 0.05, mu=1e-3)$P; toc()
  res = baseline_spectral(P, numClust=numClust)
  result_k[,i] = res
  ri_k[i] = rand.index(truth$V1, res)
}



