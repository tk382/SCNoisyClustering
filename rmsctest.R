library(matrixStats)
library(irlba)
library(SIMLR)
library(ggplot2)
library(KRLS)
library(dplyr)
library(reshape)
library(caret)
library(fossil)
data("BuettnerFlorian")
X = BuettnerFlorian$in_X
truth = BuettnerFlorian$true_labs
numClust = BuettnerFlorian$n_clust
dim(X)
reshuffle = sample(1:nrow(X))
X1 = X[reshuffle[1:2996],]
X2 = X[reshuffle[2997:5992],]
X3 = X[reshuffle[5993:8989],]

#randomly split X's features to 3 pieces
T1 = constructKernel(t(X1), 5)
T2 = constructKernel(t(X2), 5)
T3 = constructKernel(t(X3), 5)
diag(T3) = colMeans(T3)
TT = melt(T3)
ggplot(TT, aes(x=X1, y=X2, fill=value)) + geom_tile()+scale_fill_gradient(high='indianred', low='antiquewhite1')
T = array(0,dim=c(182,182,3))
T[,,1] = T1; T[,,2] = T2; T[,,3] = T3

lambda = 0.005
verbose=TRUE
eps = 1e-9
max_iter = 300
rm(T1, T2, T3)


P = RMSC(T, 0.05, verbose=TRUE)
ggplot(melt(P), aes(x=X1, y=X2, fill=value))+geom_tile()+scale_fill_gradient(high='indianred', low='antiquewhite1', lim=c(0.00, 0.03))
res = baseline_spectral(P, numClust=3)
plot(truth$V1 ~ res)
table(truth$V1, res)


rand.index(truth$V1, res)


cbind(truth$V1, res)

for (i in 1:length(res)){
  if(res[i]==3){res[i] <- 2}
  else if(res[i]==1){res[i] <-3}
  else{res[i] <- 1}
}


table(res, truth$V1)
