#' For large version, save the large kernel matrices in the hard drive
#'
#' @param X Data matrix to construct kernel from
#' @param dir Directory to save the matrix, default as the current directory
#' @param kernel_type Type of distance measure to use : euclidean, pearson, spearman, or combined
#' @param klist Kernel parameters k
#' @param sigmalist Kernel parameters sigma
#' @param verbose Show progress if set to TRUE
#'
#' @export
save.sparse.kernel = function(X,
                            dir = ".",
                            kernel_type = kernel_type,
                            klist = seq(15,25,by=5),
                            sigmalist=seq(1,2,by=0.2),
                            verbose = FALSE){
  ii = 1
  filenames = c()
  if(kernel_type %in% c("pearson", "combined")){
    ker = "pearson"
    if(verbose){print("building pearson correlation matrix..")}
    diff1 = 1-cor(as.matrix(X), method = "pearson")
    if(verbose){print("      now saving following kernels")}
    for (kk in klist){
      for (ss in sigmalist){
        P = get_kernel_matrix(X, diff1, kk, ss)
        filename = paste0(dir,'/',ker, ii, "_SLSL_kernels_pearson.mtx")
        filenames = c(filenames, filename)
        writeMM(P, filename, col.names=FALSE, row.names=FALSE, quote = FALSE)
        if(verbose){print(paste0("           ", filename))}
        ii = ii + 1
      }
    }
    rm(diff1); gc()
  }
  ii=1
  if(kernel_type %in% c("euclidean", "combined")){
    ker = "euclidean"
    if(verbose){print("building euclidean distance matrix..")}
    diff2 = as.matrix(dist(t(X)))
    if(verbose){print("      now saving following kernels")}
    for (kk in klist){
      for (ss in sigmalist){
        P = get_kernel_matrix(X, diff2, kk, ss)
        filename = paste0(dir,'/',ker, ii, "_SLSL_kernels.mtx")
        filenames = c(filenames, filename)
        # write.table(as.matrix(P), filename, col.names=FALSE, row.names=FALSE, quote = FALSE)
        writeMM(P, file=filename)
        if(verbose){print(paste0("           ", filename))}
        ii = ii+1
      }
    }
    rm(diff2); gc()
  }
  ii=1
  if(kernel_type %in% c("spearman", "combined")){
    ker = "spearman"
    if(verbose){print("building euclidean distance matrix..")}
    diff3 = 1-cor(as.matrix(X), method = "spearman")
    if(verbose){print("      now saving following kernels")}
    for (kk in klist){
      for (ss in sigmalist){
        P = get_kernel_matrix(X, diff3, kk, ss)
        P[is.na(P)] = 0
        filename = paste0(dir,'/',ker, ii, "_SLSL_kernels.txt")
        filenames = c(filenames, filename)
        writeMM(P, file=filename)
        if(verbose){print(paste0("           ", filename))}
        ii = ii+1
      }
    }
    rm(diff3); gc()
  }

  return(filenames)
}
