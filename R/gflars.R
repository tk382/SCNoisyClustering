###################################
# Y: observed matrix (p \times n)
# theta: unobserved matrix (p \times n) of true values
# phi: unknown vector of p probe effects
# delta1 = theta[1, ]
# delta2[i, ] = (theta[i + 1, ]c - theta[i, ]) / wts[i] for i = 1 : (p - 1)

gflars <- function(Y, wts = defaultWeights(nrow(Y)), steps = min(dim(Y)) - 1) {

  p <- dim(Y)[1]; n <- dim(Y)[2]

  ##wts <- rep(1, p - 1)
  ##wts <- sqrt(p / (1 : (p - 1)) / ((p - 1) : 1))

  ############################################################
  phi <- rep(1, p); phi.old <- phi + 1e+10; err.phi <- 1e+10
  loop <- 0; theta1 <- NULL
  while ((loop < 20) & (err.phi > 1e-5)) {

    theta1 <- t(phi / p) %*% Y; dim(theta1) <- NULL

    phi <- Y %*% (theta1 / sum(theta1^2))
    phi <- phi / sqrt(mean(phi^2)); dim(phi) <- NULL

    err.phi <- sqrt(sum((phi.old - phi)^2) / p)
    if (err.phi > 1e-5) {
      phi.old <- phi
      rss.0 <- sum((Y - phi %*% t(theta1))^2) / n / p / 2
      loop <- loop + 1
    }
    #print(c(loop, err.phi))
  }

  df.0 <- n + p - 1

  output <- list()
  output[[1]] <- list(theta = rep(1, p) %x% t(theta1), phi = phi, rss = rss.0,
                      loop = loop, df.aic = df.0, df.bic.y = df.0, df.bic.k = df.0, change.p = numeric(0))
  #print(c("k = ", 1))
  #########################################
  for (k in 1 : (steps - 1)) {
    ##phi <- rep(1, p)
    phi.old <- phi + 1e+10; err.phi <- 1e+10; loop <- 0

    while ((loop < 10) & (err.phi > 1e-5)) {

      grp.reg <- doGFLars(Y, k, phi = phi)

      bkp <- grp.reg$bkp
      order.bkp <- order(bkp)
      bkp <- bkp[order.bkp]
      change.p <- bkp + 1

      delta2 <- matrix(0, p - 1, n)
      delta2[bkp, ] <- grp.reg$value[[k]][order.bkp, ]

      temp <- cumsum(phi[p : 1]^2)[p : 1][-1]
      delta2.up <- delta2 * wts
      delta1 <- t(phi) %*% Y / p - t(temp) %*% delta2.up / p

      delta <- rbind(delta1, delta2.up)
      theta <- delta
      for (i in 2 : p) theta[i, ] <- theta[i - 1, ] + delta[i, ]

      change.p.set <- list()
      count <- 1; change.p.set[[1]] <- 1
      for (i in 2 : p) {
        if (is.element(i, change.p)) {
          count <- count + 1
          change.p.set[[count]] <- i
        }
        if (!is.element(i, change.p)) change.p.set[[count]] <- c(change.p.set[[count]], i)
      }

      for (g in 1 : length(change.p.set)) {
        cp.set.g <- change.p.set[[g]]
        l.g <- length(cp.set.g)
        #if (l.g == 1) phi[cp.set.g] <- 1 tae update: needs to preserve the sign
        if (l.g == 1) phi[cp.set.g] <- sign(phi[cp.set.g])
        if (l.g > 1) {
          theta.g <- theta[cp.set.g[1], ]
          phi[cp.set.g] <- Y[cp.set.g, ] %*% (theta.g / sum(theta.g^2))
          phi[cp.set.g] <- phi[cp.set.g] / sqrt(mean(phi[cp.set.g]^2))
        }
      }

      df.aic   <- length(change.p.set) * n + p - length(change.p.set) + 3 * length(change.p)
      df.bic.y <- length(change.p.set) * n + p - length(change.p.set) + length(change.p)
      df.bic.k <- length(change.p.set) * n + p - length(change.p.set) + 2 * length(change.p)

      rss <- sum((Y - theta * phi)^2) / n / p / 2

      err.phi <- sqrt(sum((phi.old - phi)^2) / p)
      print(err.phi)
      if (err.phi > 1e-5) {
        phi.old <- phi
        loop <- loop + 1
      }
        output[[k + 1]] <- list(theta = theta, phi = phi, rss = rss, loop = loop,
                                df.aic = df.aic, df.bic.y = df.bic.y, df.bic.k = df.bic.k, change.p = change.p)
#      }  #tae edit : output is updated when error is small
      #print(c(err.phi, loop))
    }
    #print(c("k = ", k + 1))
  }

  rss.path <- df.path.aic <- df.path.bic.y <- df.path.bic.k <- NULL
  for (k in 1 : length(output)) {
    rss.path <- c(rss.path, output[[k]]$rss)
    df.path.aic   <- c(df.path.aic, output[[k]]$df.aic)
    df.path.bic.y <- c(df.path.bic.y, output[[k]]$df.bic.y)
    df.path.bic.k <- c(df.path.bic.k, output[[k]]$df.bic.k)
  }
  ##rss.path[df.path == (p * n + 2 * (p - 1))] <- max(rss.path)

  #  BIC.path.aic   <- log(rss.path) + 2 * df.path.aic / p
  #  BIC.path.bic.y <- log(rss.path) + log(p) * df.path.bic.y / p
  #  BIC.path.bic.k <- log(rss.path) + log(p) * df.path.bic.k / p

  BIC.path.aic   <- log(rss.path) + 2 * df.path.aic / p / n
  #BIC.path.bic.y <- log(rss.path) + log(p * n) * df.path.bic.y / p / n
  #BIC.path.bic.k <- log(rss.path) + log(p * n) * df.path.bic.k / p / n

  #should be log(p) instead fo log(p*n)
  BIC.path.bic.y <- log(rss.path) + log(p) * df.path.bic.y / p / n
  BIC.path.bic.k <- log(rss.path) + log(p) * df.path.bic.k / p / n

  ans.aic   <- output[[which.min(BIC.path.aic)]]
  ans.bic.y <- output[[which.min(BIC.path.bic.y)]]
  ans.bic.k <- output[[which.min(BIC.path.bic.k)]]

  return(list(ans.aic = ans.aic, ans.bic.y = ans.bic.y, ans.bic.k = ans.bic.k, output = output))
}
