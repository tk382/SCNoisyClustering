heat2 = function(S, pheno){
  library(reshape2)
  rownames(S) = paste0('R',1:nrow(S))
  colnames(S) = paste0('C', 1:ncol(S))
  mx = melt(S)
  colnames(mx) = c('Var1','Var2','value')
  mx$Var1 = factor(mx$Var1, levels = rownames(S))
  mx$Var2 = factor(mx$Var2, levels = colnames(S))
  ggplot(mx, aes(x=Var1, y=Var2,fill=value)) +
    geom_tile()+
    scale_fill_gradient(low="azure2", high="black")
}
