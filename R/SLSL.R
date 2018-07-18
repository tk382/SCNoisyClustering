SLSL = function(X,
                numClust = NA,
                ref = NA,
                core = NA,
                k = NA,
                log = T,
                filter = F,
                filter_p1 = 0.9,
                filter_p2 = 0,
                correct_detection_rate = F,
                kernel_type = "combined",
                klist = seq(15,25,by=5),
                sigmalist=seq(1,2,by=0.2),
                tau = 5,
                gamma = 0,
                verbose=FALSE,
                measuretime = FALSE,
                warning=TRUE
                ){
  # X           : each column is one cell, each row is one gene
  # k           : number of neighbors for kernel computation and network diffusion
  # numClust    : number of clusters
  # kernel_type : possible options: "pearson", "euclidean", "spearman", "combined"# klist       : kernel parameters
  # sigmalist   : kernel parameters


  if(ncol(X)>1200 & warning & is.na(ref)){
    stop("We detected more than 1,200 cells, and system might crash due to memory requirement.
         We recommend LSLSL function for large matrices.
         If you'd like to use SLSL anyway, set warning=FALSE")
  }

  if(length(ref) > 1){

    out = SLSL_ref(X = X,
                   ref = ref,
                   core = core,
                   log = log,
                   numClust = numClust)
    return(out)
  }else{
    if(is.na(k)){
      k = max(10,ncol(X)/20)
    }

    if(!is.matrix(X)){
      X = as.matrix(X)
    }

    # start measuring time
    t1 = Sys.time()

    # Filter the genes
    if(filter){
      X = genefilter(X, filter_p1, filter_p2)
    }

    t2 = Sys.time()    #t2 - t1 : gene filtering

    if(log){X = log(X+1)}

    if(correct_detection_rate & kernel_type %in% c('euclidean','combined')){
      det = colSums(X!=0) / nrow(X)
      det2 = qr(det)
      X = t(qr.resid(det2, t(X)))
    }

    X = scale(X)

    # Construct kernel..
    if(verbose){print('constructing kernel..')}
    if(kernel_type=="pearson"){
      diff = 1-cor(as.matrix(X), method = "pearson")
      P = corr_kernel_c(t(X),klist,sigmalist,k)
    }else if(kernel_type=="euclidean"){
      P = dist_kernel_c(t(X), klist, sigmalist, k)
    }else if(kernel_type=="spearman"){
      diff = 1-cor(as.matrix(X), method = "spearman")
      P    = rank_kernel_c(t(X), diff, klist, sigmalist, k)
    }else if(kernel_type=="combined"){
      diff1 = 1-cor(as.matrix(X), method = "pearson")
      diff3 = 1-cor(as.matrix(X), method = "spearman")
      P1 = corr_kernel_c(t(X), diff1, klist, sigmalist, k)
      P2 = dist_kernel_c(t(X), klist, sigmalist, k)
      P3 = rank_kernel_c(t(X), diff3, klist, sigmalist, k)
      P  = abind(P1, P2, P3)
      rm(P1, P2, P3)
    }else{
      stop("kernel_type must be one of 'pearson','euclidean','spearman', or,'combined'")
    }

    t3 = Sys.time()

    if(verbose){
      print('optimizing..')
    }

    res = sparse_scaledlasso_c(P, tau, gamma, verbose=verbose)
    sigmas = res$sigma

    t4 = Sys.time()

    if(verbose){print('network diffusion..')}

    S = network.diffusion(res$S, k)

    t5 = Sys.time()

    if(verbose){print('dimension reduction..')}

    if(is.na(numClust)){
      if(verbose){print('determining cluster number..')}
      numClust= getClustNum(S)
    }

    tmp = tsne_spectral(S, numClust)
    t6 = Sys.time()
    if(measuretime){cat(paste('gene filter & numClust',
                              as.difftime(round(t2-t1,2), units="secs"),
                              'seconds \nconstructing kernel',
                              as.difftime(round(t3-t2,2), units="secs"),
                              'seconds \nget S',
                              as.difftime(round(t4-t3,2), units="secs"),
                              'seconds \nnetwork diffusion',
                              as.difftime(round(t5-t4,2), units="secs"),
                              'seconds\ntsne',
                              as.difftime(round(t6-t5,2), units="secs"),
                              'seconds'))}

    return(list(S=S, result = tmp$cluster, tsne=tmp$tsne, sigma = sigmas))
  }

}
