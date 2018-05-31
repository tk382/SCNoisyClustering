#ensemble type method for large scale clustering


ssl_setup()
#
load('data/large_Zeisel.RData')
X = Zelsel$in_X

#load('data/1_Kolod.RData')
#X = Test_2_Kolod$in_X
nn = ncol(X)
shuffle = sample(1:ncol(X))
X = X[, shuffle]
lab = Zelsel$true_labs[shuffle,1]
#lab = Test_2_Kolod$true_labs[shuffle,1]
k = 9
rm(Zelsel)
#rm(Test_2_Kolod)
X2 = list()
lablist = list()
division = 10
skip = round(ncol(X) / division)
indvector = 1:skip; indlist = list()
for (i in 1:(division-1)){
  indlist[[i]] = (1+skip*(i-1)):(skip*i)
  lablist[[i]] = lab[indvector]
  X2[[i]] = X[,indvector]
  X = X[,-indvector]
  lab = lab[-indvector]
}
X2[[division]] = X; lablist[[division]] = lab
indlist[[division]] = (indlist[[division-1]][length(indlist[[division-1]])]+1):nn
rm(X); rm(lab)
gc()

groupslist = combn(1:length(X2), 2)
N = ncol(groupslist)
cluster_result = matrix(NA, nn, N)


for (n in 1:N){
  t = Sys.time()
  groups = groupslist[,n]
  tmpX = cbind(X2[[groups[1]]], X2[[groups[2]]])
  #zeros = apply(tmpX, 1, function(x) sum(x==0))
  #remove = which(zeros > ncol(tmpX)*0.95)
  t = Sys.time()
  estimates = ssl_wrapper(tmpX, numClust = 9, verbose=FALSE, measuretime=F)
  print(Sys.time()-t) #13 seconds in one round
  cluster_result[c(indlist[[groups[1]]], indlist[[groups[2]]]),n] = estimates$result
  print(Sys.time()-t)
  gc()
  print(groups)
}


final = array(0, dim=c(nrow(cluster_result), ncol(cluster_result),1,1))
final[,,1,1]= cluster_result
k=9
dimnames(final)[[1]] = paste0('R',1:nrow(cluster_result))
dimnames(final)[[2]] = paste0('C',1:ncol(cluster_result))
dimnames(final)[[3]] = paste0("k")
dimnames(final)[[4]] = as.character(k)

final2 = CSPA(final, 9)
adj.rand.index(final2, lab)


ensemble_example_result = list(shuffle = shuffle,
              X = X2,
              lablist = lablist,
              result = cluster_result)
save(ensemble_example_result, file='ensemble_example_result.Rdata')


load('ensemble_example_result.Rdata')
cluster_result = ensemble_example_result$result
final = array(0, dim=c(nrow(cluster_result), ncol(cluster_result),1,1))
final[,,1,1]= cluster_result
k=9
dimnames(final)[[1]] = paste0('R',1:nrow(cluster_result))
dimnames(final)[[2]] = paste0('C',1:ncol(cluster_result))
dimnames(final)[[3]] = paste0("k")
dimnames(final)[[4]] = as.character(k)

final2 = CSPA(final, 9)
compare(vec, final2, method = "nmi" )
