#' Wrapper function for similarity learning algorithm. For details, refer to the manuscript.
#' We recommend users to start with default setting for filter, kernel_type, klist, sigmalist, tau, gamma, and k.
#' The only critical parameter is log. If X has been log-transformed and log is not set to FALSE, it will produce an error.
#'
#'
#' @param X Expression level matrix in count, cells in columns, genes in rows
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
LSLSL = function(X,
                 seed = NA,
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
                 sigmalist=seq(1,2,by=0.5),
                 tau = 5,
                 gamma = 0,
                 verbose=F){
  #CSPA:  Cluster-based Similarity Partitioning Algorithm
  #LCE:   Linkage Clustering Ensemble
  #majority_voting :  majority voting

  if(is.na(seed)){set.seed(sample(1:1000))}

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

  #shuffle the data (so that labels are mixed up)
  if(shuffle){
    shuffle = sample(1:ncol(X)); X = X[, shuffle]
  }else{
    shuffle = 1:ncol(X)
  }

  filenames = savesparsekernel(X = X,
                               kernel_type = kernel_type,
                               klist = klist,
                               sigmalist = sigmalist,
                               verbose = verbose)

  #split X into multiple data sets
  nn = ncol(X)
  division = round(nn/1000)
  skip = floor(nn/division)
  indvector = 1:skip; indlist = list()
  for (i in 1:(division-1)){
    indlist[[i]] = (1+skip*(i-1)):(skip*i)
  }
  indlist[[division]] = (indlist[[division-1]][length(indlist[[division-1]])]+1):nn

  #form the pairs
  gl = combn(1:length(indlist), 2)
  N = ncol(gl)
  if(is.na(core)){core = detectCores()-1}
  if (core < 1) {
    core = 1
  }
  if(core>1){
    if(verbose){print("working on the following subsets:")}
    myfun = function(l){
      groups = gl[,l]
      if(verbose){print(paste0("     ", groups))}
      tmpind = c(indlist[[groups[1]]],
                 indlist[[groups[2]]])
      P = rep(list(Matrix(0, nrow=length(tmpind), ncol = length(tmpind), sparse=T)),
              length(filenames))
      for (i in length(filenames)){
        # tmp = as.matrix(setDF(fread(filenames[i], select=tmpind))[tmpind, ])
        tmp = readMM(filenames[i])[tmpind,tmpind]
        # tmp2 = Matrix(tmp, sparse = TRUE)
        # P[[i]] = tmp2
        P[[i]] = tmp
      }
      res = sparse_scaledlasso_list_c(P=P,
                                      n=length(P),
                                      tau = tau,
                                      gamma = gamma,
                                      verbose = FALSE)
      sigmas = res$sigma
      S = network.diffusion(res$S, k)
      if(is.na(numClust)){
        numClust= getClustNum(S)
      }
      tmp = tsne_spectral(S, numClust)
      return(list(result = tmp$cluster, groups = groups))
    }
    cl = makeCluster(core, type = "FORK", outfile = "")
    res = parLapply(cl, 1:N, myfun)
    final = matrix(NA, nn, N)
    for (i in 1:N){
      item = res[[i]]
      group = item$groups
      inds = c(indlist[[group[1]]], indlist[[group[2]]])
      final[shuffle[inds], i] = item$result
    }
    stopCluster(cl)
  }else{
    if(verbose){print("working on the following subsets:")}
    #don't use parallelization
    myfun = function(l){
      groups = gl[,l]
      if(verbose){print(paste0("        ", groups))}
      tmpind = c(indlist[[groups[1]]],
                 indlist[[groups[2]]])
      P = rep(list(Matrix(0, nrow=length(tmpind), ncol = length(tmpind), sparse=T)),
              length(filenames))
      for (i in length(filenames)){
        tmp = as.matrix(setDF(fread(filenames[i], select=tmpind))[tmpind, ])
        tmp2 = Matrix(tmp, sparse = TRUE)
        P[[i]] = tmp2
      }
      res = sparse_scaledlasso_list_c(P=P,
                                      n=length(P),
                                      tau = tau,
                                      gamma = gamma,
                                      verbose = FALSE)
      sigmas = res$sigma
      S = network.diffusion(res$S, k)
      if(is.na(numClust)){
        numClust= getClustNum(S)
      }
      tmp = tsne_spectral(S, numClust)
      return(list(result = tmp$cluster, groups = groups))
    }
    res = list()
    for (l in 1:N){
      res[[l]] = myfun(l)
    }
    if(verbose){print("saving result..")}
    final = matrix(NA, nn, N)
    for (i in 1:N){
      item = res[[i]]
      group = item$groups
      inds = c(indlist[[group[1]]], indlist[[group[2]]])
      final[shuffle[inds], i] = item$result
    }
  }

  k = numClust
  # if(is.na(k)){
  #   k = length(unique(as.numeric(final)))
  # }
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
    E = impute_missing(out, t(orig), k)
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
  return(list(array = final, result = result))
}
