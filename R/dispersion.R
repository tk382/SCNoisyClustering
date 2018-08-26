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
dispersion = function(X, bins=20){

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
  # median_disp     = df %>% group_by(quartile) %>% summarize(mediandisp = median(disp))
  # newdf           = full_join(df, median_disp, by='quartile')
  # newdf           = newdf %>% mutate(num = disp - mediandisp)
  # tmp             = newdf %>% group_by(quartile) %>%
  #   summarize(denom=median(abs(disp-mediandisp)))
  # newdf           = full_join(newdf, tmp, by='quartile')
  # newdf           = newdf %>% mutate(normalized_disp = num/denom)
  # newdf$genenames = genenames
  return(newdf)
}
