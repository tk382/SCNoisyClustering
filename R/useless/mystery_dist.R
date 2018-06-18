# compute and returns the multiple kernel
mystery_dist = function(x, k = 10, allk_input=seq(10,30,by=5), sigma_input=seq(2,1,by=-0.1)){
  N = dim(x)[1] #number of cells
  KK = 0
  sigma = sigma_input

  # compute and sort Diff
  Diff = dist2(x)^2
  Diff_sort = t(apply(Diff,MARGIN=2,FUN=sort))

  # compute the combined kernels
  allk = allk_input
  D_Kernels = list()
  for (l in 1:length(allk)){
    if(allk[l] < nrow(x)-1){
      TT  = apply(Diff_sort[, 2:(allk[l]+1)], 1, mean)
      TT  = matrix(TT, nrow=length(TT), ncol=1)
      Sig = apply(array(0,c(nrow(TT),ncol(TT))),MARGIN=1,FUN=function(x){x=TT[,1]})
      Sig = (Sig + t(Sig))/2
      Sig_valid = array(0, c(nrow(Sig), ncol(Sig)))
      Sig_valid[which(Sig > .Machine$double.eps, arr.ind=TRUE)] = 1
      Sig = Sig * Sig_valid + .Machine$double.eps
      for (j in 1:length(sigma)){
        KK = KK+1
        W = dnorm(Diff, 0, sigma[j]*Sig)
        D_Kernels[[KK]] = (W + t(W))/2
      }
    }
  }

  for (i in 1:length(D_Kernels)){
    K = D_Kernels[[i]]
    dinv = 1/sqrt(diag(K)+1)
    G = K * (dinv %*% t(dinv))
    G1 = apply(array(0,c(length(diag(G)),length(diag(G)))),MARGIN=2,FUN=function(x) {x=diag(G)})
    G2 = t(G1)
    #G1 and G2 are diagonal elements of G repeated in each row and column
    D_Kernels_tmp = (G1 + G2 - 2*G)/2 #<- magic step : discovers 3rd cluster! what is this doing??
    D_Kernels_tmp2 = D_Kernels_tmp - diag(diag(D_Kernels_tmp))
    D_Kernels[[i]] = D_Kernels_tmp
  }

  return(D_Kernels)

}
