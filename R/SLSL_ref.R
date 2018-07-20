SLSL_ref = function(X,
                    ref,
                    #knn = T,
                    #knn_keep = 10,
                    numClust = NA,
                    k = NA,
                    log = T,
                    filter = F,
                    filter_p1 = 0.9,
                    filter_p2 = 0,
                    correct_detection_rate = F,
                    kernel_type = "combined",
                    klist = seq(15,25,by=5),
                    sigmalist = seq(1,2,by=0.2),
                    tau = 5,
                    gamma = 0,
                    verbose = FALSE,
                    measuretime = FALSE,
                    warning = TRUE){

  # take subset of the data and
  # the reference set with overlapping genes
  int = base::intersect(rownames(ref), rownames(X))
  if(length(int)<30){
    stop("Not enough genes overlap with the reference data set.
         Check that row names of the reference set and rownames of your data set have gene names.")
  }


  if(verbose){print('data cleaning..')}
  ind1 = match(int, rownames(X))
  X = X[ind1, ]
  ind2 = match(int, rownames(ref))
  ref = ref[ind2, ]


  if(log){
    X = log(X+1)
  }

  if(verbose){'correcting..'}
  if(correct_detection_rate){
    det = colSums(X!=0) / nrow(X)
    det2 = qr(det)
    X = t(qr.resid(det2, t(X)))
    X = scale(X, center=T, scale=T)
  }

  #X = scale(X)

  if(verbose){print('computing projection..')}

  projection = proj_c(as.matrix(X), as.matrix(ref))

  # if(knn){
  #   ranks = apply(projection, 2, function(x) match((1:nrow(projection)),order(x, decreasing = TRUE)))
  #   projection[which(ranks>knn_keep, arr.ind=TRUE)] = 0
  #
  #   #remove cell typess that were not selected in knn
  #   zeros = rowSums(projection==0)
  #   ind = which(zeros==ncol(projection))
  #   if(length(ind)>0){projection = projection[-ind, ]}
  #   projection = scale(projection^4)
  #   colnames(projection) = colnames(X)
  #   rownames(projection) = colnames(ref)[-ind]
  # }else{
  #   #normalize
  #   projection = scale(projection^4)
  #   colnames(projection) = colnames(X)
  #   rownames(projection) = colnames(ref)
  # }



  if(is.na(numClust)){
    numClust = getClustNum(projection)
  }

  # if(verbose){print('hierarchical clustering..')}
  # result = cutree(hclust(dist(t(projection))), numClust)

  if(ncol(projection) < 5000){
    print('ref with smaller version..')
    res = SLSL2(as.matrix(projection),
               numClust = numClust,
               k = NA,
               log = F,
               filter = F,
               filter_p1 = 1,
               filter_p2 = 0,
               correct_detection_rate = F,
               kernel_type = kernel_type,
               klist = klist,
               sigmalist = sigmalist,
               tau = tau,
               gamma = gamma,
               verbose = verbose,
               measuretime = measuretime,
               warning = warning
    )
  }else{
    print('ref with large version..')
    res = LSLSL(as.matrix(projection),
                numClust = numClust,
                core = core,
                shuffle=TRUE,
                cluster_method = "CSPA",
                k = NA,
                log = F,
                filter = T,
                filter_p1 = 0.9,
                filter_p2 = 0,
                correct_detection_rate = F,
                kernel_type = "combined",
                klist = seq(15,25,by=5),
                sigmalist=seq(1,2,by=0.2),
                tau = 5,
                gamma = 0,
                verbose = verbose
                )
  }
  out = list(result = res, projection = projection)
  #out = list(SLSL_output = res, projection = projection, numGenes=length(int))
  return(out)
}
