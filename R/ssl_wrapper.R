ssl_wrapper = function(X,
                       numClust = NA,
                       k=NA,
                       klist = seq(10,30,by=5),
                       sigmalist=seq(2,1,by=-0.1),
                       filter = TRUE,
                       verbose=TRUE,
                       measuretime = FALSE
                       ){
  # X : each column is one cell, each row is one gene
  # k : number of neighbors for kernel computation and network diffusion
  # numClust : number of clusters
  # klist : kernel parameters
  # sigmalist : kernel parameters
  #
  t1 = Sys.time()
  X = genefilter(X)
  if(is.na(numClust)){
    numClust= getClustNum(X)
  }
  t2 = Sys.time()
  if(is.na(k)){k = min(10,ncol(X)/20)}
  if(verbose){print('constructing kernel..')}
  P = corr_kernel_c(t(X),klist,sigmalist,k)
  t3 = Sys.time()
  if(verbose){print('optimizing..')}
  res = sparse_scaledlasso_c(P, 5, 0, verbose=FALSE)
  t4 = Sys.time()
  if(verbose){print('network diffusion..')}
  S = network.diffusion(res$S, k)
  t5 = Sys.time()
  if(verbose){print('dimension reduction..')}
  tmp = tsne_spectral(S, numClust)
  t6 = Sys.time()
  if(measuretime){print(paste('gene filter & numClust', t2-t1,
                              'construcg kernel', t3-t2,
                              'get S', t4-t3,
                              'network diffusion', t5-t4,
                              'tsne', t6-t5))}
  return(list(S=S, result = tmp$cluster, tsne=tmp$tsne))
}
