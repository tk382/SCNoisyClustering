#' Remove genes (rows) with little information.
#' Genes that are not expressed in any of the samples, or most of the samples, and in some cases, genes that are expressed in most the samples, will be removed.
#'
#' @param X data matrix with zeros
#' @param bins number of bins to divide genes by their mean expression level

#' @return data frame with mean and variance of expression level for each gene, bins they belong to, and normalized dispersion level
#' @export
plot.dispersion = function(X,
                           genenames,
                           bins=NA, median = FALSE,
                           outliers.mean.thresh = c(30,Inf),
                           outliers.vmr.thresh = c(3,Inf)){
  if(sum(X<0) > 0){
    warning("The input matrix X should be count without log-transformation")
  }
  if(is.na(bins)){
    bins=10
  }
  genemeans  = Matrix::rowMeans(X)
  genevars   = apply(X, 1, var)
  disp       = genevars / genemeans
  quartile   = ntile(genemeans, bins)
  df         = data.frame( genemeans  = genemeans,
                           genevars   = genevars,
                           disp       = disp,
                           quartile   = quartile )
  if(!median){
    newdf = df %>%
      dplyr::group_by(quartile) %>%
      dplyr::mutate(normalized_dispersion = scale(disp))
    newdf$normalized_dispersion[!is.finite(newdf$normalized_dispersion)] = 0
    newdf$normalized_dispersion[is.na(newdf$normalized_dispersion)] = 0
  }else{
    newdf = df %>% dplyr::group_by(quartile) %>%
      dplyr::mutate(bin_median = median(disp)) %>%
      dplyr::mutate(bin_mad = mad(disp))
    newdf = newdf %>% dplyr::group_by(quartile) %>%
      dplyr::mutate(normalized_dispersion = abs(disp - bin_median) / bin_mad)
    newdf$normalized_dispersion[!is.finite(newdf$normalized_dispersion)] = 0
    newdf$normalized_dispersion[is.na(newdf$normalized_dispersion)] = 0
  }
  newdf$vmr = log(newdf$normalized_dispersion+1)
  ind = 1:nrow(df)
  outliers = which((newdf$vmr>outliers.vmr.thresh[1] &
                     newdf$vmr < outliers.vmr.thresh[2]) |
                     (newdf$genemeans > outliers.mean.thresh[1] &
                     newdf$genemeans < outliers.mean.thresh[2]))
  plot(newdf$vmr[-outliers] ~ newdf$genemeans[-outliers],
       ylab="log(vmr)",
       xlab="",
       main="log of normalized dispersion",
       ylim = c(min(newdf$vmr), max(newdf$vmr)),
       xlim = c(min(newdf$genemeans), max(newdf$genemeans)),
       cex = 0.5)
  text(newdf$vmr[outliers] ~ newdf$genemeans[outliers],
       labels = genenames[outliers],
       cex=0.7)
  return(newdf)
}

