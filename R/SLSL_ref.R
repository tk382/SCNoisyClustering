#' Wrapper function for similarity learning algorithm using reference panel. For details, refer to the manuscript.
#' We recommend users to start with default setting for filter, kernel_type, klist, sigmalist, tau, gamma, and k.
#' The only critical parameter is log. If X has been log-transformed and log is not set to FALSE, it will produce an error.
#'
#'
#' @param X Expression level matrix in count, cells in columns, genes in rows. The gene names should be the row names of X.
#' @param ref Reference panel, cells in columns, genes in rows. The gene names in the row names.
#' @param numClust Number of Clusters. If unknown, NA. If set to NA, it will be estimated with minimum 3.
#' @param log If set to TRUE, SLSL will take log from X. If it is already log-transformed, set FALSE.
#' @param filter If set to TRUE, remove rows with little information
#' @param filter_p1 Upper threshold for percentage of zeros in each row
#' @param filter_p2 Lower threshold for percentage of zeros in each row
#' @param correct_detection_rate If set to TRUE, detection rate will be regressed out
#' @param kernel_type distance measure to use: one of "euclidean", "pearson", "spearman", or "combined"
#' @param klist Kernel parameters setting.
#' @param sigmalist Kernel parameters setting.
#' @param tau Regularization parameter for L1 norm
#' @param gamma Regularization parameter for Frobenius norm
#' @param k Number of neighbors to use for knn procedure.
#' @param verbose To show progress, set to TRUE
#' @param measuretime Measure the time it takes for each step
#' @param warning Warns if the data size is too large. To ignore, set to FALSE.
#' @examples
#'
#' sysdata = SCNoisyClustering::reference_panel
#' ref = reference_panel$cell

#' load('../data/unnecessary_in_building/7_Chu_celltype.RData')
#' X = Chu_celltype$X
#' chu_ref = SLSL_ref(X = as.matrix(X),
#'                   ref = ref,
#'                   numClust = 7,
#'                   verbose = T)
#'
#'
SLSL_ref = function(X,
                    ref,
                    numClust = NA,
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
                    k = NA,
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

  # if(is.na(knn_keep)){
  #   knn_keep = round(ncol(ref)/3)
  # }


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
  }else{
    X = scale(X)
  }



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
    res = SLSL(as.matrix(projection),
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
  return(out)
}
