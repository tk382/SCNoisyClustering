library(matrixStats)
library(irlba)
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
library(microbenchmark)
##MUST LOAD BEFORE COMPILING RCPP##
library(inline)
library(quadprog)
##MUST LOAD##
source("../SIMLR/R/SIMLR.R")
source("../SIMLR/R/compute.multiple.kernel.R")
source("../SIMLR/R/network.diffusion.R")
source("../SIMLR/R/utils.simlr.R")
source("../SIMLR/R/tsne.R")
dyn.load("../SIMLR/R/projsplx_R.so")


set.seed(1)
R.utils::sourceDirectory("R/", modifiedOnly=FALSE)
Rcpp::sourceCpp('src/SCNoisyClustering.cpp')

#first data
data_index = 1; data_name = 'mECS'
load('../data/Test_1_mECS.RData')
dat = Test_1_mECS
rm(Test_1_mECS); gc()
X = dat$in_X
truth = dat$true_labs
numClust = dat$n_clust

#test all methods
set.seed(1)
allk = seq(2,10,by=2); sigma=seq(2,0.5,by=-0.5)
P = multiple_kernel_tae(t(X), allk, sigma)

#res = RMSC(P = P, tau = 1e-1, gamma = 0, lambda = rep(1/20, 20), mu=1, verbose=TRUE)
res = nuclear_objective_c(P, tau=1, mu=1e-3, verbose=TRUE)
ggplot(melt(res$S), aes(x=X1, y=X2, fill=value)) + geom_tile() + scale_color_gradient()


K = makeK(X, allk = seq(2,20,by=2), sigma = 1, bestind=1)
set.seed(1)
w1 = objective_1_c(P, numClust, mu=1e-3, rho = 1.9, max_iter=100, eps=1e-7, verbose=TRUE)
makeplotandsave(w1$S, "mECS_w")
set.seed(1)
g1 = objective_2_c(P, numClust, mu=1e-3, rho=1.9, max_iter=100, eps=1e-7, verbose=TRUE)
makeplotandsave(g1$S, "mECS_g")
set.seed(1)
m1 = RMSC_c(K, lambda=0.005, verbose=TRUE)
#makeplotandsave(g1$S, "mECS_m")
set.seed(1)
s1 = SIMLR(X, numClust)
#makeplotandsave(s1$S, "mECS_s")
set.seed(1)
k1 = kmeans(t(X), numClust)
rm(P, K); gc()

#see results
#change the first line and assign the result of interest in res_1
res_1 = res #w1, g1, m1, s1, k1
S_1 = res_1$S; #f_1 = res_1$f; w_1 = res_1$w
clust_1 = baseline_spectral(S_1, numClust=numClust)
ri_1 = rand.index(truth$V1, clust_1)
nmi_1 = compare(truth$V1, clust_1, method='nmi')
print(c(ri_1, nmi_1))
table(truth$V1, clust_1)
ggplot(melt(S_1), aes(x=X1, y=X2, fill=value)) + geom_tile() + scale_color_gradient()
clust_p1 = k1$cluster
ri_p1 = rand.index(truth$V1, clust_p1)
nmi_p1 = compare(truth$V1, clust_p1, method="nmi")
print(c(ri_p1, nmi_p1))

#pollen data
data_name = 'Pollen'
load('../data/Test_3_Pollen.RData')
dat = Test_3_Pollen
rm(Test_3_Pollen); gc();
X = dat$in_X
truth = dat$true_labs
numClust = dat$n_clust
ind = sort(truth$V1, index.return=TRUE)$ix

#test all methods
#test all methods
set.seed(1)
allk = seq(2,20,by=2); sigma=seq(2,0.2,by=-0.2)
P = multiple_kernel_tae(t(X), allk, sigma)
set.seed(1)
w2 = objective_1_c(P, numClust, mu=1e-3, rho = 1.9, max_iter=100, eps=1e-4, verbose=TRUE)
makeplotandsave(w2$S, "Pollen_w")
set.seed(1)
#w22 = objective_1_c_one_aux(P, numClust, mu=1e-3, rho=1.9, max_iter = 100, eps = 1e-4, verbose=TRUE)

g2 = objective_2_c(P, numClust, mu=1e-3, rho=1.9, max_iter=100, eps=1e-4, verbose=TRUE)
makeplotandsave(g2$S, "Pollen_g")
set.seed(1)
m2 = RMSC_c(K, lambda=0.005, verbose=TRUE)
#makeplotandsave(g1$S, "mECS_m")
set.seed(1)
s2 = SIMLR(X, numClust, allk, sigma)
#makeplotandsave(s1$S, "mECS_s")
set.seed(1)
#kmeans directly on data
k2 = kmeans(t(X), numClust, nstart = 20)
rm(P, K); gc()
#pca on data and keep nClust eigenvectors
p2 = irlba(X, numClust)
tmpdf_p = data.frame(truth = as.factor(truth$V1), v1 = p2$v[,1], v2 = p2$v[,2], v3 = p2$v[,3], v4 = p2$v[,4])
scatterplot3d(x = tmpdf_p$v1, y = tmpdf_p$v2, z = tmpdf_p$v3, color=tmpdf_p$truth)

ww2 = irlba(s2$S, numClust)
tmpdf_w = data.frame(truth = as.factor(truth$V1), v1 = ww2$v[,1], v2 = ww2$v[,2], v3 = ww2$v[,3], v4 = ww2$v[,4])
scatterplot3d(x = tmpdf_w$v2, y=tmpdf_w$v3, z=tmpdf_w$v4, color=tmpdf_w$truth)

res_2 = g2 #w2, g2, m2, s2, k2
S_2 = res_2$S
f_2 = res_2$f
w_2 = res_2$w
clust_2 = baseline_spectral(S_2, numClust=numClust)
ri_2 = adj.rand.index(truth$V1, clust_2)
nmi_2 = compare(truth$V1, clust_2, method='nmi')
print(c(ri_2, nmi_2))
#table(truth$V1, clust_2)
ggplot(melt(S_2[ind,ind]), aes(x=X1, y=X2, fill=value)) + geom_tile() + scale_fill_gradient(high="indianred",low="skyblue")

clust_p2 = k2$cluster
ri_p2 = rand.index(truth$V1, clust_p2)
nmi_p2 = compare(truth$V1, clust_p2, method="nmi")
print(c(ri_p2, nmi_p2))

#inspect allk and sigma
df = data.frame(k = rep(allk, each=10), sigma=rep(sigma, 10))
df$weight = g2$sigma
ggplot(df, aes(x=k, y=sigma, fill=weight)) + geom_tile() + scale_fill_gradient()




par(mfrow=c(1,2))
g1 = ggplot(melt(P[ind,ind,5]), aes(x=X1, y=X2, fill=value)) + geom_tile() + scale_color_gradient()
g2 = ggplot(melt(ddi[ind,ind]), aes(x=X1, y=X2, fill=value))+geom_tile() + scale_color_gradient()
grid.arrange(g1, g2, nrow=1)
K = constructKernel(t(X), 10)
ggplot(melt(K), aes(x=X1, y=X2, fill=value))+geom_tile() + scale_color_gradient()


#Chu et al. cell types
X = read.table('../data/GSE75748_sc_cell_type_ec_cleaned.txt', header=TRUE, stringsAsFactors = FALSE)
truth = data.frame(V1 = as.factor(sapply(strsplit(colnames(X), '_'), function(x) x[1])))
truth$V1 =as.integer(truth$V1)
numClust = 7
gc()

#test all methods
#test all methods
set.seed(1)
allk = seq(2,20,by=2); sigma=seq(1, 0.1, by= -0.1)
P = multiple_kernel_tae(t(X), allk, sigma)
#K = makeK(X, allk = seq(2,20,by=2), sigma = 1, bestind=1)
set.seed(1)
w3 = objective_1_c(P, numClust, mu=1e-3, rho = 1.9, max_iter=100, eps=1e-7, verbose=TRUE)
#makeplotandsave(w1$S, "mECS_w")
set.seed(1)
g3 = objective_2_c(P, numClust, mu=1e-3, rho=1.9, max_iter=100, eps=1e-7, verbose=TRUE)
#makeplotandsave(g1$S, "mECS_g")
set.seed(1)
m3 = RMSC_c(K, lambda=0.005, verbose=TRUE)
#makeplotandsave(g1$S, "mECS_m")
set.seed(1)
s3 = SIMLR(X, numClust, allk, sigma)
#makeplotandsave(s1$S, "mECS_s")
set.seed(1)
k3 = kmeans(t(X), numClust)
rm(P, K); gc()

res_3 = w3 #w3, g3, m3, s3, k3
S_3 = res_3$S
f_3 = res_3$f
w_3 = res_3$w
clust_3 = baseline_spectral(S_3, numClust=numClust)
ri_3 = adj.rand.index(truth$V1, clust_3)
nmi_3 = compare(truth$V1, clust_3, method='nmi')
print(c(ri_3, nmi_3))
table(truth$V1, clust_3)
ggplot(melt((S_3+t(S_3))/2), aes(x=X1, y=X2, fill=value)) + geom_tile() + scale_color_gradient()
clust_p3 = k3$cluster
ri_p3 = rand.index(truth$V1, clust_p3)
nmi_p3 = compare(truth$V1, clust_p3, method="nmi")
print(c(ri_p3, nmi_p3))





#Chu et al. 2nd experiment
X = read.table('cleaned_celltype_exp2.txt', header=TRUE, stringsAsFactors = FALSE)
truth = data.frame(V1 = as.factor(sapply(strsplit(colnames(X), '_'), function(x) x[1])))
truth$V1 =as.integer(truth$V1)
numClust = length(unique(truth$V1))

#test all methods
set.seed(1)
allk = seq(10,30,by=2); sigma=seq(2,1,by=-0.25)
#P = makeP(X, allk, sigma)
P = multiple_kernel_tae(t(X), allk, sigma)
K = makeK(X, allk, sigma=1, bestind=1)
set.seed(1)
w4 = objective_1_c(P, numClust, eps=1e-9, verbose=TRUE, max_iter=100)
makeplotandsave(w4$S, "Chu2_w")
set.seed(1)
g4 = objective_2_c(P, numClust, eps=1e-9, verbose=TRUE, max_iter=100)
makeplotandsave(g4$S, "Chu2_g")
set.seed(1)
m4 = RMSC_c(K, lambda=0.005, verbose=TRUE)
#makeplotandsave(g4$S, "Chu2_m")
set.seed(1)
s4 = SIMLR(X, numClust)
makeplotandsave(s4$S, "Chu2_s")
set.seed(1)
k4 = kmeans(t(X), numClust)
gc()

res_4 = w4    #w4, g4, m4, s4, k4
S_4 = res_4$S
f_4 = res_4$f
w_4 = res_4$w
clust_4 = baseline_spectral(S_4, numClust=numClust)
ri_4 = rand.index(truth$V1, clust_4)
nmi_4 = compare(truth$V1, clust_4, method='nmi')
print(c(ri_4, nmi_4))
table(truth$V1, clust_4)
ggplot(melt(S_4), aes(x=X1, y=X2, fill=value)) + geom_tile() + scale_color_gradient()
clust_p4 = k4$cluster
ri_p4 = rand.index(truth$V1, clust_p4)
nmi_p4 = compare(truth$V1, clust_p4, method="nmi")
print(c(ri_p4, nmi_p4))



# microbenchmark(objective_2(as.matrix(X), numClust, makeP(X), eps=1e-6, max_iter=100),
#                objective_2_c(makeP(X), numClust, eps=1e-6, max_iter=100),
#                times=10)
#
# microbenchmark(objective_1(as.matrix(X), numClust, makeP(X), eps=1e-6, max_iter=10),
#                objective_1_c(makeP(X), numClust, eps=1e-6, max_iter=10),
#                times = 10)
#


# P = makeP(X)
# microbenchmark(objective_1_c(P, numClust, eps=1e-6, max_iter=100),
#                objective_2_c(P, numClust, eps=1e-6, max_iter=100),
#                times=10L)
#objective 2 is faster if we wait until convergence, probably due to quadprog
