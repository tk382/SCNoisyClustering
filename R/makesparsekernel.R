#' Estimate the optimal number of clusters using eigengap
#'
#' @param S estimated similarity matrix.
#' @examples
#' #create positive definite symmetric matrix
#' X = matrix(rnorm(50), nrow = 10)
#' S= t(X) %*% X
#' getClustNum(S)
#'
makesparsekernel = function(X,
                            kernel_type = kernel_type,
                            klist = seq(15,25,by=5),
                            sigmalist=seq(1,2,by=0.2),
                            verbose = FALSE){
  ii = 1
  if(kernel_type%in% c("pearson", "combined")){
    P1 = rep(list(Matrix(0, nrow=ncol(X), ncol = ncol(X), sparse=T)),
             length(klist) * length(sigmalist))
    diff1 = 1-cor(as.matrix(X), method = "pearson")
    ii=1
    for (kk in klist){
      for (ss in sigmalist){
        P1[[ii]] = get_kernel_matrix(X, diff1, kk, ss)
        ii = ii+1
      }
    }
    rm(diff1); gc()
  }
  ii=1
  if(kernel_type %in% c("euclidean", "combined")){
    P2 = rep(list(Matrix(0, nrow=ncol(X), ncol = ncol(X), sparse=T)),
             length(klist) * length(sigmalist))
    diff2 = as.matrix(dist(t(X)))
    #diff2 = dist_c(t(X))
    for (kk in klist){
      for (ss in sigmalist){
        P2[[ii]] = get_kernel_matrix(X, diff2, kk, ss)
        ii = ii+1
      }
    }
    rm(diff2); gc()
  }
  ii=1
  if(kernel_type %in% c("spearman", "combined")){
    P3 = rep(list(Matrix(0, nrow=ncol(X), ncol = ncol(X), sparse=T)),
             length(klist) * length(sigmalist))
    diff3 = 1-cor(as.matrix(X), method = "spearman")
    for (kk in klist){
      for (ss in sigmalist){
        tmp = get_kernel_matrix(X, diff3, kk, ss)
        tmp[is.na(tmp)] = 0
        P3[[ii]] = tmp
        ii = ii+1
      }
    }
    rm(diff3); gc()
  }

  if(kernel_type == "pearson"){
    return(P1)
  }else if(kernel_type=="euclidean"){
    return(P2)
  }else if(kernel_type=="spearman"){
    return(P3)
  }else{
    return(c(P1,P2,P3))
  }
}
