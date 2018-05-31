large_ssl = function(X, numClust, shuffle = TRUE, core = NA){
  library(parallel)
  nn = ncol(X)
  if(shuffle){
    shuffle = sample(1:ncol(X))
    X = X[, shuffle]
  }
  X2 = list()
  division = round(nn/500)
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

  gl = combn(1:length(X2), 2)
  N = ncol(gl)
  cluster_result = matrix(NA, nn, N)

  if(is.na(core)){core = detectCores()-1}

  if (core < 1 || is.na(core) || is.null(core)) {
    core = 1
  }
  cl = makeCluster(core)

  clusterEvalQ(cl, {library(Matrix)})
  res = list()
  res = parLapply(cl, 1:N,
                         fun=function(l,
                                      groupslist_fun = gl,
                                      X_fun=X2,
                                      numClust_fun = numClust){
    print(class(groupslist_fun))
    groups = groupslist_fun[,l]
    print(groups)
    tmpX = cbind(X_fun[[groups[1]]], X_fun[[groups[2]]])
    estimates = ssl_wrapper(tmpX, numClust_fun, verbose=TRUE, measuretime=F)
    return(estimates$result)
  }
  )
  stopCluster(cl)


}
