dist2 = function( x, c = NA ) {

  # set the parameters for x
  if(is.na(c)) {
    c = x
  }

  # compute the dimension
  n1 = nrow(x)
  d1 = ncol(x)
  n2 = nrow(c)
  d2 = ncol(c)
  if(d1!=d2) {
    stop("Data dimension does not match dimension of centres.")
  }

  dist = t(rep(1,n2) %*% t(rowSums(x^2))) +
    rep(1,n1) %*% t(rowSums(c^2)) - 2 * (x %*% t(c))

  return(dist)

}
