#' Show distribution of normalized dispersion for each gene
#'
#'
#' @param X data matrix
#' @param genenames name of the genes
#' @param bins number of bins to divide genes by their mean expression level
#' @param median whether to use median instead of mean, MAD instead of SD
#' @param outliers.mean.thresh x-axis threshold to show outliergene names instead of simple plotting
#' @param outliers.vmr.thresh y-axis threshold to show outlier gene names instead of simple plotting
#' @export
plot_dispersion = function(X,
                           genenames,
                           bins=10, median = TRUE,
                           outliers.mean.thresh = c(30,Inf),
                           outliers.vmr.thresh = c(5,Inf)){
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
  if(length(outliers)>0){
    plot(newdf$vmr[-outliers] ~ newdf$genemeans[-outliers],
       ylab="log(vmr)",
       xlab="",
       main="log of normalized dispersion",
       ylim = c(min(newdf$vmr, na.rm=T), max(newdf$vmr, na.rm=T)),
       xlim = c(min(newdf$genemeans, na.rm=T), max(newdf$genemeans, na.rm=T)),
       cex = 0.5)
    text(newdf$vmr[outliers] ~ newdf$genemeans[outliers],
         labels = genenames[outliers],
         cex=0.7)
  }else{
    plot(newdf$vmr ~ newdf$genemeans,
         ylab="log(vmr)", xlab="", main="log of normalized dispersion")
  }

  return(newdf)
}

