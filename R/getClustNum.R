getClustNum = function(X){
  s = svd(X)
  eigengap = abs(diff(s$d))
  c = max(3, which.max(eigengap[-1])+1)
  return(c)
}
