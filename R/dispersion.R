#' Remove genes (rows) with little information.
#' Genes that are not expressed in any of the samples, or most of the samples, and in some cases, genes that are expressed in most the samples, will be removed.
#'
#' @param X data matrix with zeros
#' @param bins number of bins to divide genes by their mean expression level

#' @return data frame with mean and variance of expression level for each gene, bins they belong to, and normalized dispersion level
#' @examples
#' X = matrix(sample(0:10, size = 10000, replace = TRUE, prob = c(0.9, rep(0.1/10, 10))), nrow = 200) #create expression level matrix
#' newX = genefilter(X, 0.9, 0)
#' dim(newX)
#' @export
dispersion = function(X, bins=20, robust = FALSE){
  if(!robust){
  if(sum(X<0) > 0){
    warning("The input matrix X should be count without log-transformation")
  }

  genemeans  = Matrix::rowMeans(X)
  genevars   = apply(X, 1, var)
  disp       = genevars / genemeans
  quartile   = ntile(genemeans, bins)
  df         = data.frame( genemeans  = genemeans,
                           genevars   = genevars,
                           disp       = disp,
                           quartile   = quartile )
  newdf = df %>% group_by(quartile) %>% mutate(z = scale(disp))
  return(newdf)
  }else{
    df = data.frame(mean=Matrix::colMeans(X),
                   cv=apply(X,2,sd)/Matrix::colMeans(X),
                   var=apply(X,2,var))
    df$dispersion = with(df,var/mean)
    df$mean_bin = with(df,cut(mean,breaks=c(-Inf,quantile(mean,seq(0.1,1,0.05)),Inf)))
    var_by_bin = ddply(df,"mean_bin",function(x) {
      data.frame(bin_median = median(x$dispersion),
                 bin_mad = mad(x$dispersion))
    })
    df$bin_disp_median = var_by_bin$bin_median[match(df$mean_bin,var_by_bin$mean_bin)]
    df$bin_disp_mad = var_by_bin$bin_mad[match(df$mean_bin,var_by_bin$mean_bin)]
    df$dispersion_norm = with(df,abs(dispersion-bin_disp_median)/bin_disp_mad)
    return(df)
  }



}
