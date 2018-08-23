#' Remove genes (rows) with little information.
#' Genes that are not expressed in any of the samples, or most of the samples, and in some cases, genes that are expressed in most the samples, will be removed.
#'
#'
#' @param X data matrix with zeros
#' @param p1 upper threshold of the percentage of zeros in each row
#' @param p2 lower threshold of the percentage of zeros in each row
#' @return The data matrix with rows with too many or too few zeros removed
#' @examples
#' X = matrix(sample(0:10, size = 10000, replace = TRUE, prob = c(0.9, rep(0.1/10, 10))), nrow = 200) #create expression level matrix
#' newX = genefilter(X, 0.9, 0)
#' dim(newX)
dispersion = function(X){
  # zeros = rowSums(X==0)
  # remove = which(zeros > (p1 * ncol(X)))
  # remove2 = c(remove, which(zeros < (p2 * ncol(X))))
  # if(length(remove2)>0){X = X[-remove2,]}
  # return(X)
  if(sum(X<0) > 0){
    warning("The input matrix X should be count without log-transformation")
  }
  genemeans  = rowMeans(X)
  genevars   = apply(X, 1, var)
  disp       = genevars / genemeans
  quartile   = ntile(genemeans, 20)
  df         = data.frame( genemeans  = genemeans,
                           genevars   = genevars,
                           disp       = disp,
                           quartile   = quartile )
  median_disp = df %>% group_by(quartile) %>% summarize(mediandisp = median(disp))
  newdf       = full_join(df, median_disp, by='quartile')
  newdf       = newdf %>% mutate(num = disp - mediandisp)
  tmp         = newdf %>% group_by(quartile) %>%
    summarize(denom=median(abs(disp-mediandisp)))
  newdf       = full_join(newdf, tmp, by='quartile')
  newdf       = newdf %>% mutate(normalized_disp = num/denom)
  rownames(newdf) = rownames(X)
  return(newdf[,c(1,2,3,4,8)])
}
