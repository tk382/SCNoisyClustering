LSLSL = function(X,
                 numClust = NA,
                 kernel_type = "combined",
                 core = NA,
                 shuffle=TRUE,
                 verbose=F){

  library(parallel)

  nn = ncol(X)

  #shuffle the data (so that labels are mixed up)
  if(shuffle){
    shuffle = sample(1:ncol(X))
    X = X[, shuffle]
  }

  #split X into multiple data sets
  X2 = list()
  division = round(nn/400)
  skip = floor(nn/division)
  indvector = 1:skip; indlist = list()
  for (i in 1:(division-1)){
    indlist[[i]] = (1+skip*(i-1)):(skip*i)
    X2[[i]] = X[,indvector]
    X = X[,-indvector]
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
    cl = makeCluster(core, type="FORK")
    myfun = function(l){
      groups = gl[,l]
      if(verbose){print(groups)}
      tmpX = cbind(X2[[groups[1]]], X2[[groups[2]]])
      tmpX = log(tmpX+1)
      estimates = SLSL(X = tmpX, numClust = numClust, kernel_type = kernel_type, verbose=T)
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

  result = CSPA(out, k)

  return(result)
}
