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
  if(!robust){
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
  return(newdf)


}
