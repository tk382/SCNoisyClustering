SLSL = function(X,
                numClust = NA,
                ref = "none",
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
                verbose=FALSE,
                measuretime = FALSE
                ){
  # X           : each column is one cell, each row is one gene
  # k           : number of neighbors for kernel computation and network diffusion
  # numClust    : number of clusters
  # kernel_type : possible options: "pearson", "euclidean", "spearman", "combined"# klist       : kernel parameters
  # sigmalist   : kernel parameters

  if(ref %in% c("tissue", "cell")){
    if(ref=="tissue"){
      load('data/sysdata.rda')
      ref = sysdata$GlobalPanel[[2]]
    }else if(ref=="cell"){
      load('data/sysdata.rda')
      ref = sysdata$GlobalPanel[[1]]
    }
    X = log(X+1)
    out = SLSL_ref(X, ref, numClust)
    return(out)

  }else if(ref=="none"){
    if(is.na(k)){
      k = max(10,ncol(X)/20)
    }

    if(!is.matrix(X)){
      X = as.matrix(X)
    }

    # start measuring time
    t1 = Sys.time()

    # Filter the genes
    X = genefilter(X, filter_p1, filter_p2)

    t2 = Sys.time()    #t2 - t1 : gene filtering

    if(log){X = log(X+1)}

    X2 = X

    if(correct_detection_rate & kernel_type %in% c('euclidean','combined')){
      zeros = apply(X2, 2, function(x) sum(x==0))
      pc1 = irlba(X2, 1)$v[,1]
      mod = lm(pc1~zeros)
      if(tidy(mod)[2,5] < 0.05){
        if(verbose){print('correcting for the detection rate..')}
        scalar = sqrt(ncol(X)-1)
        X3 = scale(t(X2)) / scalar
        x = scale(zeros) / scalar
        # x = scale(pc1)/scalar
        X2 = t(X3 - x  %*% (t(x) %*% X3)) * scalar
      }
      X2 = scale(X2, scale=F, center=T)
    }

    # Construct kernel..
    if(verbose){print('constructing kernel..')}
    if(kernel_type=="pearson"){
      P = corr_kernel_c(t(X),klist,sigmalist,k)
    }else if(kernel_type=="euclidean"){
      P = dist_kernel_c(t(X2), klist, sigmalist, k)
    }else if(kernel_type=="spearman"){
      diff = 1-cor(as.matrix(X), method = "spearman")
      P    = rank_kernel_c(t(X), diff, klist, sigmalist, k)
    }else if(kernel_type=="combined"){
      diff = 1-cor(as.matrix(X), method = "spearman")
      P1 = corr_kernel_c(t(X), klist, sigmalist, k)
      P2 = dist_kernel_c(t(X2), klist, sigmalist, k)
      P3 = rank_kernel_c(t(X), diff, klist, sigmalist, k)
      P  = abind(P1, P2, P3)
      rm(P1, P2, P3)
    }else{
      stop("kernel_type must be one of 'pearson','euclidean','spearman', or,'combined'")
    }

    t3 = Sys.time()

    if(verbose){print('optimizing..')}

    res = sparse_scaledlasso_c(P, tau, gamma, verbose=F)
    sigmas = res$sigma

    t4 = Sys.time()

    if(verbose){print('network diffusion..')}

    S = network.diffusion(res$S, k)

    t5 = Sys.time()

    if(verbose){print('dimension reduction..')}

    if(is.na(numClust)){
      print('determining cluster number..')
      numClust= getClustNum(S)
    }

    tmp = tsne_spectral(S, numClust)
    t6 = Sys.time()
    if(measuretime){cat(paste('gene filter & numClust',
                              as.difftime(round(t2-t1,2), units="secs"),
                              '\nconstructing kernel',
                              as.difftime(round(t3-t2,2), units="secs"),
                              '\nget S',
                              as.difftime(round(t4-t3,2), units="secs"),
                              '\nnetwork diffusion',
                              as.difftime(round(t5-t4,2), units="secs"),
                              '\ntsne',
                              as.difftime(round(t6-t5,2), units="secs")))}

    return(list(S=S, result = tmp$cluster, tsne=tmp$tsne))
  }else{
    stop("ref should be one of 'none', 'tissue', or 'cell'.")
  }

}
