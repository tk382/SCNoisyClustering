makegif = function(S, thresh, name){
  library(caTools)
  giffile = array(0, dim=c(ncol(S),ncol(S),length(thresh)))
  i = 1
  for (th in thresh){
    S1 = as.matrix(makeA_thresh(S, th))
    giffile[,,i] = S1
    i = i+1
  }
  write.gif(giffile, name, scale = "always")
}
