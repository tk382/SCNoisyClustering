makeA_thresh = function(S, lambda){
  S[S>=lambda] = 1
  S[S<1] = 0
  return(S)
}
