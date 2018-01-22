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

#which data set to test?

data_index = 1; data_name = 'mECS'
#data_index = 3; data_name = 'Pollen'
#data_index = 4; data_name = 'Usoskin'

#load data and define X, true labels, and number of clusters
datname = paste0('Test_',data_index, '_',data_name)
load(paste0('../data/',datname,'.Rdata'))
assign("dat", get(datname))
X = dat$in_X
truth = dat$true_labs
numClust = dat$n_clust
rm(data_index, data_name, datname)
rm(Test_1_mECS)

#Chu et al. cell types
# X = read.table('../data/GSE75748_sc_cell_type_ec_cleaned.txt', header=TRUE, stringsAsFactors = FALSE)
# truth = data.frame(V1 = as.factor(sapply(strsplit(colnames(X), '_'), function(x) x[1])))
# truth$V1 =as.integer(truth$V1)
# numClust = 7

res_1 = objective_1(X, numClust, verbose=TRUE)
S_1 = res_1$S
f_1 = res_1$f
w_1 = res_1$w
clust_1 = baseline_spectral(S_1, numClust=numClust)
ri_1 = rand.index(truth$V1, clust_1)
nmi_1 = compare(truth$V1, clust_1, method='nmi')
table(truth$V1, clust_1)
ggplot(melt(S_1), aes(x=X1, y=X2, fill=value)) + geom_tile() + scale_color_gradient()


res_2 = objective_2(X, numClust, verbose=TRUE)
S_2 = res_2$S
f_2 = res_2$f
sigma_2 = res_2$sigma
clust_2 = baseline_spectral(S_2, numClust=numClust)
ri_2 = rand.index(truth$V1, clust_2)
nmi_2 = compare(truth$V1, clust_2, method='nmi')
table(truth$V1, clust_2)
ggplot(melt(S_2), aes(x=X1, y=X2, fill=value)) + geom_tile() + scale_color_gradient()


res_mv = RMSC_multiview(X, 0.005, 3, verbose=TRUE)
S_mv = res_mv$P
clust_mv = baseline_spectral(S_mv, numClust = numClust)
ri_mv = rand.index(truth$V1, clust_mv)
nmi_mv = compare(truth$V1, clust_mv, method='nmi')
ggplot(melt(S_mv), aes(x=X1, y=X2, fill=value)) + geom_tile() + scale_color_gradient()


source("../SIMLR/R/SIMLR.R")
source("../SIMLR/R/compute.multiple.kernel.R")
source("../SIMLR/R/network.diffusion.R")
source("../SIMLR/R/utils.simlr.R")
source("../SIMLR/R/tsne.R")
dyn.load("../SIMLR/R/projsplx_R.so")
res_s = SIMLR(X, 3)
S_s = res_s$S
w_s = res_s$alphaK
clust_s = baseline_spectral(S_s, numClust = numClust)
ri_s = rand.index(truth$V1, clust_s)
nmi_s = compare(truth$V1, clust_s, method = 'nmi')
ggplot(melt(S_s), aes(x=X1, y=X2, fill=value)) + geom_tile() + scale_color_gradient()

gc()

S_1_data1  = ggplot(melt(S_1), aes(x=X1, y=X2, fill=value)) + geom_tile() + scale_fill_gradient(high='white',low='indianred4')
S_2_data1  = ggplot(melt(S_2), aes(x=X1, y=X2, fill=value)) + geom_tile() + scale_fill_gradient(high='white',low='indianred4')
S_mv_data1 = ggplot(melt(S_mv), aes(x=X1, y=X2, fill=value)) + geom_tile() + scale_fill_gradient(high='white',low='indianred4')
S_s_data1  = ggplot(melt(S_s), aes(x=X1, y=X2, fill=value)) + geom_tile() + scale_fill_gradient(high='white',low='indianred4')


library(gridExtra)
combined = grid.arrange(S_1_data1,S_2_data1,S_mv_data1,S_s_data1,
             S_1_data2,S_2_data2,S_mv_data2,S_s_data2,
             S_1_data3,S_2_data3,S_mv_data3,S_s_data3,
             nrow=3)


ggsave(file='dumb.png', combined, width=20)
