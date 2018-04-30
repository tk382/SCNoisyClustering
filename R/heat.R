heat = function(X){
  ggplot(melt(X), aes(x=X1, y=X2,fill=value)) +
    geom_tile()+
    scale_fill_gradient(low="azure2", high="black")
}
