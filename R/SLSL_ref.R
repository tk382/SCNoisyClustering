SLSL_ref = function(X,
                    ref,
                    numClust = NA,
                    k = NA,
                    log = T,
                    filter = F,
                    filter_p1 = 0.9,
                    filter_p2 = 0,
                    correct_detection_rate = F,
                    kernel_type = "combined",
                    klist = seq(15,25,by=5),
                    sigmalist = seq(1,2,by=0.2),
                    tau = 5,
                    gamma = 0,
                    verbose = FALSE,
                    measuretime = FALSE,
                    warning = TRUE){

  # take subset of the data and
  # the reference set with overlapping genes
  int = base::intersect(rownames(ref), rownames(X))
  if(length(int)<100){
    stop("Not enough genes overlap with the reference data set.
         Check that row names of the reference set and rownames of your data set have gene names.")
  }

  ind1 = match(int, rownames(X))
  X = X[ind1, ]
  ind2 = match(int, rownames(ref))
  ref = ref[ind2, ]


  if(log){
    X = log(X+1)
  }

  X2 = X
  if(correct_detection_rate){
    det = colSums(X!=0) / nrow(X)
    det2 = qr(det)
    X2 = t(qr.resid(det2, t(X)))
    X2 = scale(X2, center=T, scale=T)
  }else{
    X2 = X
  }

  projection = proj_c(as.matrix(X2), as.matrix(ref))
  #t(t(X) %*% as.matrix(ref))/nrow(X)

  projection = scale(projection^4)

  res = SLSL(as.matrix(projection),
             numClust = NA,
             ref = NA,
             k = NA,
             log = F,
             filter = F,
             filter_p1 = 1,
             filter_p2 = 0,
             correct_detection_rate = F,
             kernel_type = kernel_type,
             klist = klist,
             sigmalist = sigmalist,
             tau = tau,
             gamma = gamma,
             verbose = verbose,
             measuretime = measuretime,
             warning = warning
             )
  out = list(SLSL_output = res, projection = projection, numGenes=length(int))
  return(out)
}
