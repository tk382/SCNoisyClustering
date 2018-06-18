ssl = function(X, nClust, k=10,
               filter = FALSE,
               d = 0.04,
               logtransform = FALSE,
               allk_list=seq(10,30, by=2),
               sigma_list=seq(2,1, by=-0.2)){
  #if X matrix is count, logtransform should be set to TRUE
  #X should be p by n matrix, where n is the number of cells and p is the number of genes
  #nClust argument should be removed and determing C should be implemented
  if(logtransform){X = log(X+1)}
  n = ncol(X); p = nrow(X)
  #gene filtering
  if(filter){
    print('gene filtering..')
    zerocounts = apply(X, 1, function(x) sum(x==0))
    ind = which(zerocounts < (d * n) | zerocounts > ((1-d)*n))
    X = X[-ind, ]
  }

  #construct kernel
  print('constructing kernel..')
  P = mystery_kernel(t(X), k=10, allk_list, sigma_list)

  #run the algorithm
  print('run optimization..')
  ssl = sparse_scaledlasso(P,5,0,nClust,verbose=FALSE)
  S = ssl$S

  #network diffusion
  print('network diffusion..')
  S2 = network.diffusion(S, nClust)

  #dimension reduction and clustering
  print('tsne..')
  res = tsne_spectral(S2, nClust) #can fix numEigen argument
  return(res)
}
