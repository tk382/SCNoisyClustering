#' Wrapper function for similarity learning algorithm. For details, refer to the manuscript.
#' We recommend users to start with default setting for filter, kernel_type, klist, sigmalist, tau, gamma, and k.
#' The only critical parameter is log. If X has been log-transformed and log is not set to FALSE, it will produce an error.
#'
#'
#' @param X Expression level matrix in count, cells in columns, genes in rows
#' @param niter Number of randomly selected subset to run unsupervised SLSL
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
#' #create sample count matrix
#' X = matrix(rpois(100000, 1), nrow = 1000)
#' result = SLSL(X)
#'
#' @export
LSLSL_supervised = function(X,
                            seed = NA,
                            niter = 5,
                            numClust = NA,
                            core = NA,
                            shuffle=TRUE,
                            cluster_method = "CSPA",
                            k=NA,
                            log = T,
                            filter = T,
                            filter_p1 = 0.9,
                            filter_p2 = 0,
                            correct_detection_rate = F,
                            kernel_type = "combined",
                            klist = seq(15,25,by=5),
                            sigmalist=seq(1,2,by=0.2),
                            tau = 5,
                            gamma = 0,
                            verbose=F){
  #CSPA:  Cluster-based Similarity Partitioning Algorithm
  #LCE:   Linkage Clustering Ensemble
  #majority_voting :  majority voting

  if(is.na(seed)){set.seed(sample(1:1000))}else{set.seed(seed)}
  if(!kernel_type %in% c("pearson","euclidean","spearman","combined")){
    stop("kernel_type must be one of 'pearson','euclidean','spearman', or,'combined'")
  }
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
    X = scale(X)
  }

  X = t(scale(t(X)))
  final = matrix(NA, ncol(X), niter)
  selected = matrix(0, 2000, niter)
  for (i in 1:niter){
    select = sample(1:ncol(X), size = 2000)
    selected[,i] = select
    tmpX = X[, select]
    otherX = X[, -select]
    res = SLSL(tmpX, log = F, filter = F, correct_detection_rate = F,
               kernel_type = kernel_type, klist = klist, sigmalist = sigmalist,
               verbose=verbose, measuretime = FALSE, warning=FALSE)
    final[select,i] = res$result

    proj = proj_c(tmpX, otherX)
    nc = length(unique(res$result))
    ref = matrix(NA, nrow(tmpX), nc)
    for (k in 1:nc){
      ref[,k] = rowMeans(tmpX[,res$result==k])
    }
    proj = proj_c(otherX, ref)
    rest_label = apply(proj, 2, which.max)
    final[-select, i] = rest_label
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
    if(verbose){print("consensus clustering: CSPA..")}
    result = CSPA(out, k)
  }else{
    E = impute_missing(out, t(X), k)
    if(cluster_method=="LCE"){
      if(verbose){print("consensus clustering: LCE..")}
      result = LCE(E, k)
    }else if(cluster_method=="majority_voting"){
      if(verbose){print("consensus clustering: Majority Voting..")}
      result = majority_voting(E, k)
    }else{
      if(verbose){print("consensus clustering: K Modes..")}
      result = k_modes(E, k)
    }
  }
  return(list(array = final, result = result, selected = selected))
}
