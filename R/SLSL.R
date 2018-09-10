#' Wrapper function for similarity learning algorithm. For details, refer to the manuscript.
#' We recommend users to start with default setting for filter, kernel_type, klist, sigmalist, tau, gamma, and k.
#' The only critical parameter is log. If X has been log-transformed and log is not set to FALSE, it will produce an error.
#'
#'
#' @param X Expression level matrix in count, cells in columns, genes in rows
#' @param numClust Number of Clusters. If unknown, NA. If set to NA, it will be estimated with minimum 3.
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
SLSL = function(X,
                numClust = NA,
                kernel_type = "combined",
                klist = NA,
                sigmalist=seq(1,2,by=0.5),
                tau = 5,
                gamma = 0,
                k = NA,
                verbose=FALSE,
                warning=TRUE){

  if(length(klist)==1 & is.na(klist[1])){
    tmp = max(round(ncol(X)/10), 10)
    klist = c(round(tmp)*0.05, round(tmp)*0.1, round(tmp)*0.15, round(tmp)*0.2)
  }
  if(ncol(X)>5000 & warning){
    stop("We detected more than 5,000 cells, and system might crash due to memory requirement.
         We recommend LSLSL function for large matrices.
         If you'd like to use SLSL anyway with enough RAM, set warning=FALSE")
  }
  if(!kernel_type %in% c("pearson","euclidean","spearman","combined")){
    stop("kernel_type must be one of 'pearson','euclidean','spearman',or,'combined'")
  }
  if(is.na(k)){
    # k = max(10,ncol(X)/20)
    k = round(mean(klist))
  }

  if(!is.matrix(X)){
    X = as.matrix(X)
  }

  # start measuring time
  t1 = Sys.time()

  # Construct kernel..
  if(verbose){print('constructing kernel..')}
  P = makesparsekernel(X,
                       kernel_type,
                       klist,
                       sigmalist)

  t3 = Sys.time()

  if(verbose){
    print('optimizing..')
  }
  nn = length(P)
  res = sparse_scaledlasso_list_c(P = P,
                                  n = nn,
                                  tau = tau,
                                  gamma = gamma,
                                  verbose = verbose)
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
  if(verbose){print("dimension reduction..")}
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

  return(list(S=Matrix(S, sparse=TRUE),
              result = tmp$cluster,
              tsne = tmp$tsne,
              sigma = sigmas))
  }
