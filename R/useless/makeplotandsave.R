makeplotandsave = function(S, title){
  plo = ggplot(melt(S), aes(x=X1, y=X2, fill=value)) +
    geom_tile() +
    scale_fill_gradient(high='indianred', low='antiquewhite1')+
    ggtitle(title)+
    xlab("")+ylab("")+labs(fill="")
  ggsave(file=paste0(title,'.png'), plo, width=6, height=5)
}
