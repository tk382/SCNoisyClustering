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
load(paste0('../',datname,'.Rdata'))
assign("dat", get(datname))
X = dat$in_X
truth = dat$true_labs
numClust = dat$n_clust
rm(data_index, data_name, datname)
rm(Test_1_mECS)

#build feature matrices
n = 10
sigma = seq(5,10, length=n)
P = array(0, dim=c(ncol(X), ncol(X), n))
for (i in 1:n){
  Ki = constructKernel(t(X), sigma[i])
  D = diag(1/rowSums(Ki))
  P[,,i] = D %*% Ki
}
plot((sort(as.numeric(P[,,1]),decreasing=TRUE))[-(1:182)], type = 'l')
lines((sort(as.numeric(P[,,3]),decreasing=TRUE))[-(1:182)], col='red')
lines((sort(as.numeric(P[,,5]),decreasing=TRUE))[-(1:182)], col='blue')
lines((sort(as.numeric(P[,,7]),decreasing=TRUE))[-(1:182)], col='green')
lines((sort(as.numeric(P[,,9]),decreasing=TRUE))[-(1:182)], col='orange')

remove_diag_T = function(i){
  return(as.numeric(T[,,i])[-(seq(0,181)*182 + 1:182)])
}

plot(remove_diag_T(1))
points(remove_diag_T(3),col='red')


#run RMSC
res = RMSC(P, 0.1, 0, rep(1/10,10),mu=0.001, verbose=TRUE)
E = res$E
S = res$S
f = res$f
lambda = res$lambda
#make clusters based on P
res1 = baseline_spectral(S, numClust=numClust)
plot(svd(S)$v[,1:2])
table(truth$V1, res1)
rand.index(truth$V1, res1)
compare(truth$V1, res1, method='nmi')

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

ggplot(melt(S), aes(x=X1, y=X2, fill=value)) + geom_tile() +
  scale_color_gradient()



E1 = E[,,2]
diag(E1) = colMeans(E1)
ggplot(melt(E1), aes(x=X1, y=X2, fill=value)) +
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



