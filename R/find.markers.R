#' Remove genes (rows) with little information.
#' Genes that are not expressed in any of the samples, or most of the samples, and in some cases, genes that are expressed in most the samples, will be removed.
#'
#' @param X log-transformed expression level matrix with cells at columns, genes at rows
#' @param genenames gene names
#' @param cluster : clustering result
#' @param tsne : tsne result
#' @param top.n : number of top genes to return for each cluster
#' @param plot.n : number of genes to plot for each cluster
#' @return expression level matrix containing only highly variable genes
#' @export
find.markers = function(X, genenames, cluster, tsne, top.n = 100, plot.n=3){
  numClust = length(unique(cluster))
  li = list()
  plotlist = list()
  out1 = data.frame(tsne1 = tsne[,1], tsne2 = tsne[,2])
  for (i in 1:numClust){
    ind = which(cluster==i)
    p = rep(0, nrow(X))
    for (k in 1:nrow(X)){
      p[k] = suppressWarnings(wilcox.test(X[k,ind], X[k,-ind],
                                          alternative = "greater"))$p.value
    }
    p = -log10(p)
    sp = sort(p, index.return=TRUE, decreasing=TRUE)
    df = data.frame(genenames = genenames[sp$ix[1:top.n]],
                        p.values = p[sp$ix[1:top.n]])
    colnames(df) = c(paste0('clust',i,'_genenames'),
                     paste0('clust',i,'_log10p'))
    li[[i]] = df
    markers = as.character(df[1:plot.n, 1])
    ind = match(markers, genenames)
    out2 = t(X[ind, ])
    colnames(out2) = markers
    out2 = melt(out2)
    out = cbind(out2, out1)
    plotlist[[i]] = ggplot(out, aes(x=tsne1, y=tsne2, col=value)) + geom_point(size = 0.2, alpha = 0.5) +
      facet_grid(.~X2) +
      scale_color_gradient(high="black", low="antiquewhite1") +
      ggtitle(paste0('cluster', i))
  }
  return(list(plots = plotlist, markers = li))
}

# get_marker_genes <- function(dataset, labels) {
#   res <- apply(dataset, 1, get_auroc, labels = labels)
#   res <- data.frame(matrix(unlist(res), ncol = 3, byrow = T))
#   colnames(res) <- c("auroc", "clusts", "pvalue")
#   res$pvalue <- p.adjust(res$pvalue)
#   return(res)
# }
# get_auroc = function(gene, labels) {
#   score = rank(gene)
#   ms = aggregate(score ~ labels, FUN = mean)
#   posgroup <- ms[ms$score == max(ms$score), ]$labels
#   if (length(posgroup) > 1) {
#     return(c(NA, NA, NA))
#   }
#   truth <- as.numeric(labels == posgroup)
#   # Make predictions & get auc using RCOR package.
#   pred <- prediction(score, truth)
#   val <- unlist(performance(pred, "auc")@y.values)
#   pval <- suppressWarnings(wilcox.test(score[truth == 1], score[truth == 0])$p.value)
#   return(c(val, posgroup, pval))
# }
