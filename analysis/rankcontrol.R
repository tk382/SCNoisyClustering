#one last try at rank control
i=3
S = ssl[[i]]$S
sv = svd(S)
plot(sv$d)

ri = c()
klist = c(20,50,100,150,200)
for (k in klist){
  d = sv$d
  d[-(1:k)] = 0
  newS = sv$u %*% diag(d) %*% t(sv$v)
  heat(newS)
  ri = c(ri, adj.rand.index(truth[[i]], tsne_spectral(S, numClust[[i]])))
}
plot(ri)

#doesn't work... it's unpredictable. The clustering result has nothing to do with rank.
