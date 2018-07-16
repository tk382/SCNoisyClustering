LSLSL = function(X,
                 numClust = NA,
                 core = NA,
                 shuffle=TRUE,
                 cluster_method = "CSPA",
                 k=NA,
                 log = T,
                 filter = T,
                 filter_p1 = 0.9,
                 filter_p2 = 0,
                 correct_detection_rate = T,
                 kernel_type = "combined",
                 klist = seq(15,25,by=5),
                 sigmalist=seq(1,2,by=0.2),
                 tau = 5,
                 gamma = 0,
                 verbose=F){
  #CSPA:  Cluster-based Similarity Partitioning Algorithm
  #LCE:   Linkage Clustering Ensemble
  #majority_voting :  majority voting

  orig = X
  library(parallel)

  nn = ncol(X)

  #shuffle the data (so that labels are mixed up)
  if(shuffle){
    shuffle = sample(1:ncol(X)); X = X[, shuffle]
  }else{
    shuffle = 1:ncol(X)
  }

  #split X into multiple data sets
  X2 = list()
  division = round(nn/300)
  skip = floor(nn/division)
  indvector = 1:skip; indlist = list()
  for (i in 1:(division-1)){
    indlist[[i]] = (1+skip*(i-1)):(skip*i)
    X2[[i]] = X[,indvector]; X = X[,-indvector]
  }

  X2[[division]] = X;
  indlist[[division]] = (indlist[[division-1]][length(indlist[[division-1]])]+1):nn
  rm(X);

  #form the pairs
  gl = combn(1:length(X2), 2)
  N = ncol(gl)

  if(is.na(core)){core = detectCores()-1}
  if (core < 1) {
    core = 1
  }
  if(core>1){
    #set up parallelization
    cl = makeCluster(core, type = "FORK")
    myfun = function(l){
      groups = gl[,l]
      if(verbose){print(groups)}
      tmpX = cbind(X2[[groups[1]]], X2[[groups[2]]])
      tmpX = log(tmpX+1)
      estimates = SLSL(X           = tmpX,
                       numClust    = numClust,
                       kernel_type = kernel_type,
                       k           = k,
                       filter      = filter,
                       filter_p1   = filter_p1,
                       filter_p2   = filter_p2,
                       correct_detection_rate = correct_detection_rate,
                       klist       = klist,
                       sigmalist   = sigmalist,
                       tau         = tau,
                       gamma       = gamma,
                       verbose=verbose)
      return(list(result = estimates$result, groups = groups))
    }
    res = parLapply(cl, 1:N, myfun)

    final = matrix(NA, nn, N)
    for (i in 1:N){
      item = res[[i]]
      group = item$groups
      inds = c(indlist[[group[1]]], indlist[[group[2]]])
      final[shuffle[inds], i] = item$result
    }
    stopCluster(cl)
  }else{
    #don't use parallelization
    myfun = function(l){
      groups = gl[,l]
      if(verbose){print(groups)}
      tmpX = cbind(X2[[groups[1]]], X2[[groups[2]]])
      tmpX = log(tmpX+1)
      estimates = SLSL(X = tmpX, numClust = numClust, kernel_type = kernel_type, verbose=T)
      return(list(result = estimates$result, groups = groups))
    }
    res = lapply(cl, 1:N, myfun)
    final = matrix(NA, nn, N)
    for (i in 1:N){
      item = res[[i]]
      group = item$groups
      inds = c(indlist[[group[1]]], indlist[[group[2]]])
      final[shuffle[inds], i] = item$result
    }
  }

  k = numClust
  if(is.na(k)){
    k = length(unique(as.numeric(final)))
  }

  out = array(0, dim=c(nrow(final), ncol(final),1,1))
  out[,,1,1] = final
  dimnames(out)[[1]] = paste0('R',1:nrow(final))
  dimnames(out)[[2]] = paste0('C',1:ncol(final))
  dimnames(out)[[3]] = paste0("k")
  dimnames(out)[[4]] = as.character(k)

  if(cluster_method=="CSPA"){
    result = CSPA(out, k)
  }else{

    E = impute_missing(out, t(orig), k)

    if(cluster_method=="LCE"){
    result = LCE(E, k)
    }else if(cluster_method=="majority_voting"){
      result = majority_voting(E, k)
    }else{
      result = k_modes(E, k)
    }
  }


  return(list(array = final, result = result))
}
