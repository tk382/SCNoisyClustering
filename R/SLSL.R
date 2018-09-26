#' Wrapper function for similarity learning algorithm. For details, refer to the manuscript.
#' We recommend users to start with default setting for filter, kernel_type, klist, sigmalist, tau, gamma, and k.
#' The only critical parameter is log. If X has been log-transformed and log is not set to FALSE, it will produce an error.
#'
#'
#' @param X Expression level matrix in count, cells in columns, genes in rows
#' @param plot if set to TRUE, shows the plot of the dimension reduction result
#' @param numClust Number of Clusters. If unknown, NA. If set to NA, it will be estimated with minimum 3.
#' @param kernel_type distance measure to use: one of "euclidean", "pearson", "spearman", or "combined"
#' @param klist kernel parameters setting
#' @param sigmalist kernel parameters setting.
#' @param tau Regularization parameter for L1 norm
#' @param gamma Regularization parameter for Frobenius norm
#' @param verbose To show progress, set to TRUE
#' @export
SLSL = function(X,
                plot = TRUE,
                numClust = NA,
                kernel_type = "combined",
                klist = NA,
                sigmalist=c(1,2,3),
                tau = 5,
                gamma = 0,
                cluster_method = "CSPA",
                verbose=FALSE,
                warning=TRUE){

  if(length(klist)==1 & is.na(klist[1])){
    tmp = max((ncol(X)/10), 10)
    klist = c(round(tmp/2), round(tmp), round(tmp*3/2))
  }
  if(ncol(X)>5000){
    warning("We detected more than 5,000 cells, and system might crash due to memory requirement.
            We are running LSLSL instead. Using CSPA to combine results...")
    LSLSL(X,
          plot = plot,
          numClust = numClust,
          seed = NA,
          core = NA,
          shuffle=TRUE,
          cluster_method = c("CSPA", "LCE", "majority_voting"),
          kernel_type = kernel_type,
          klist = klist,
          sigmalist=sigmalist,
          tau = tau,
          gamma = gamma,
          verbose=verbose)
  }
  if(!kernel_type %in% c("pearson","euclidean","spearman","combined")){
    stop("kernel_type must be one of 'pearson','euclidean','spearman',or,'combined'")
  }
  k = round(mean(klist))

  if(!is.matrix(X)){
    X = as.matrix(X)
  }

  if(verbose){print('constructing kernel..')}
  P = make_sparse_kernel(X,
                       kernel_type,
                       klist,
                       sigmalist)
  if(verbose){print('optimizing..')}
  nn = length(P)
  res = sparse_scaledlasso_list_c(P = P,
                                  n = nn,
                                  tau = tau,
                                  gamma = gamma,
                                  verbose = verbose)
  sigmas = res$sigma
  rm(P); gc();

  if(verbose){print('network diffusion..')}
  S = network_diffusion(res$S, k)

  if(is.na(numClust)){
    numClust = get_cluster_number(S)
  }

  if(verbose){print("dimension reduction..")}
  U = dimension_reduction(S)
  # rs = rowSums(S)
  # L = (diag(rs)-S)/sqrt(rs)
  # L = t(t(L)/sqrt(rs))
  # L = (L + t(L))/2
  # eigen_L = eigen(L, symmetric = TRUE)
  # U = eigen_L$vectors
  U_index = seq(ncol(U), (ncol(U)-numClust+1))
  tmp = kmeans(U[,U_index], centers=numClust, nstart=500)
  tsne = NA
  if(plot){
    tsne = Rtsne(U[,U_index], 2, pca=FALSE, perplexity=(ncol(X)-1)/20)
    plot(tsne$Y, col=rainbow(numClust)[tmp$cluster],
         ylab="tsne2", xlab="tsne1")
  }

  return(list(S=Matrix(S, sparse=TRUE),
              result = tmp$cluster,
              tsne   = tsne))
  }
