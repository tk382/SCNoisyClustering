SLSL = function(X,
                kernel_type = "combined",
                correct_detection_rate = T,
                numClust = NA,
                k=NA,
                filter = T,
                klist = seq(15,25,by=5),
                sigmalist=seq(1.3,1.8,by=0.1),
                verbose=TRUE,
                measuretime = FALSE
                ){

  # X           : each column is one cell, each row is one gene
  # kernel_type : possible options: "pearson", "euclidean", "spearman", "combined"
  # numClust    : number of clusters
  # k           : number of neighbors for kernel computation and network diffusion
  # klist       : kernel parameters
  # sigmalist   : kernel parameters


  # if k is not specified..
  if(is.na(k)){
    k = max(10,ncol(X)/20)
  }

  # if X is not a matrix..
  if(!is.matrix(X)){
    X = as.matrix(X)
  }


  # start measuring time
  t1 = Sys.time()


  # Filter the genes
  X = genefilter(X)

  t2 = Sys.time()    #t2 - t1 : gene filtering

  X = log(X+1)
  X2 = X
  if(correct_detection_rate){
    zeros = apply(X2, 2, function(x) sum(x==0))
    pc1 = irlba(X2, 1)$v[,1]
    mod = lm(pc1~zeros)
    if(tidy(mod)[2,5] < 0.05){
      x = zeros - mean(zeros)
      X2 = t(t(scale(X2)) - x  %*% t(x) %*% t(scale(X2))/sum(x^2))
    }
  }

  # Construct kernel..
  if(verbose){print('constructing kernel..')}
  if(kernel_type=="pearson"){
    P = corr_kernel_c(t(X),klist,sigmalist,k)
  }else if(kernel_type=="euclidean"){
    P = dist_kernel_c(t(X2), klist, sigmalist, k)
  }else if(kernel_type=="spearman"){
    diff = 1-cor(as.matrix(X))
    P    = rank_kernel_c(t(X), diff, klist, sigmalist, k)
  }else if(kernel_type=="combined"){
    diff = 1-cor(as.matrix(X))
    P1 = corr_kernel_c(t(X), klist, sigmalist, k)
    P2 = dist_kernel_c(t(X2), klist, sigmalist, k)
    P3 = rank_kernel_c(t(X), diff, klist, sigmalist, k)
    P  = abind(P1, P2, P3)
  }else{
    stop("kernel_type must be one of 'pearson','euclidean','spearman', or,'combined'")
  }

  t3 = Sys.time()

  if(verbose){print('optimizing..')}

  res = sparse_scaledlasso_c(P, 5, 1, verbose=verbose)

  t4 = Sys.time()

  if(verbose){print('network diffusion..')}

  S = network.diffusion(res$S, k)

  t5 = Sys.time()

  if(verbose){print('dimension reduction..')}

  if(is.na(numClust)){numClust= getClustNum(S)}

  tmp = tsne_spectral(S, numClust)
  t6 = Sys.time()
  if(measuretime){cat(paste('gene filter & numClust', t2-t1,
                              '\n construcg kernel', t3-t2,
                              '\n get S', t4-t3,
                              '\n network diffusion', t5-t4,
                              '\n tsne', t6-t5))}

  return(list(S=S, result = tmp$cluster, tsne=tmp$tsne))
}
