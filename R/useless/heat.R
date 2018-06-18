heat = function(X, title){
  mx = melt(X)
  ggplot(mx, aes(x=X1, y=X2, fill=value)) +
    geom_tile()+
    scale_fill_gradient()

}
