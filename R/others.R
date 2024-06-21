##############################################
# Matrix form
##############################################
invert.q <- function(coef) {
  out <- 1

  if (all(abs(coef) < 1)) {
    minmod <- min(Mod(polyroot(c(1, coef))))

    if (minmod <= 1) {
      out <- 0
    }
  } else {
    out <- 0
  }

  return(out)
}

pars.mat <- function(n, parsVec, norder = 1, checkinv = TRUE) {
  Check <- 1
  if (checkinv == TRUE) {
    Check <- invert.q(parsVec)
  }
  
  if (Check == 0) {
    NULL
  } else {
    Mat <- diag(n)
    for (i in 1:norder) {
      Mat <- Mat + Diag(rep(parsVec[i], n - i), k = -i)
    }
    Mat
  }
}


#' simulate realizations using INAR(3) with zero-inflated Poisson innovation and one sustained shift
#' 
#' @param n is the length
#' @param alpha is the alpha
#' @param lambda is the mean of poisson mixture
#' @param pi is the proportion of zeros
#' @param h is the start point of shift
#' @param delta is the value of the standardized shift
#' @param burnin is the length of the burn-in period
#' @export
#' @examples
#' nsim <- 100
#' burnin <- 100
#' T <- 100
#' q <- 5
#' H <- getHMatMT(T, q)
#' Y <- arima.sim(list(ar = 0.5), n = T)
#' 
#' alpha <- c(0.03083069, 0.06242601, 0.09120189)
#' lambda <- 0.239385
#' pi <- 0.1453097
#'
#' TT <- 183
#' w <- 28
#' Y <- rzinpoisinar3(TT + w, alpha, lambda, pi, ceiling(TT / 2) + w, delta = 1, burnin = burnin)
#' 
sigma.mat <- function(n, order = c(1, 0, 0), phi.vec = 0.5, theta.vec = NULL, sigma2 = 1, burn.in = 50) {
  if (order[1] == 0) {
    phiMat <- diag(n + burn.in)
  } else {
    phiMat <- pars.mat(n + burn.in, -phi.vec, norder = order[1])
  }

  if (order[3] == 0) {
    thetaMat <- diag(n + burn.in)
  } else {
    thetaMat <- pars.mat(n + burn.in, theta.vec, norder = order[3], checkinv = FALSE)
  }

  out <- solve(phiMat) %*% thetaMat %*% t(thetaMat) %*% t(solve(phiMat)) * sigma2

  gamma0 <- out[dim(out)[1], dim(out)[2]]

  if (burn.in > 0) {
    out <- out[-c(1:burn.in), -c(1:burn.in)]
  }

  list(sigma.mat = out, sqrtsigma.mat = sqrtm(out)$B, gamma0 = gamma0)
}



mu1f <- function(delta, lambda0, pi0) {
  
  mu0 <- (1 - pi0) * lambda0
  c0 <- (1 + pi0 / (1 - pi0) * mu0)
  
  b <- - (2 * mu0 + delta ^ 2 * c0)
  c <- mu0 ^ 2
  a <- 1
  d <- b ^ 2 - 4 * a * c
  
  out <- rep(NA, 2)
  
  out[1] <- -b - sqrt(d)
  out[2] <- -b + sqrt(d)
  
  out <- out / 2 / a
  
  out[out >= mu0][1]
  
}

sigma21f <- function(delta, lambda0, pi0) {
  c0 <- (1 + pi0 * lambda0)
  mu1 <- mu1f(delta, lambda0, pi0)
  c0 * mu1
}

pi1f <- function(delta, lambda0, pi0) {
  mu1 <- mu1f(delta, lambda0, pi0)
  sigma21 <- sigma21f(delta, lambda0, pi0)
  pi1 <- 1 - mu1 ^ 2 / (sigma21 - mu1 + mu1 ^ 2)
  pi1
}

lambda1f <- function(delta, lambda0, pi0) {
  mu1 <- mu1f(delta, lambda0, pi0)
  pi1 <- pi1f(delta, lambda0, pi0)
  lambda1 <- mu1 / (1 - pi1)
  lambda1
}

#' simulate realizations using INAR(3) with zero-inflated Poisson innovation and one sustained shift
#' 
#' @param n is the length
#' @param alpha is the alpha
#' @param lambda is the mean of poisson mixture
#' @param pi is the proportion of zeros
#' @param h is the start point of shift
#' @param delta is the value of the standardized shift
#' @param burnin is the length of the burn-in period
#' @export
#' @examples
#' nsim <- 100
#' burnin <- 100
#' T <- 100
#' q <- 5
#' H <- getHMatMT(T, q)
#' Y <- arima.sim(list(ar = 0.5), n = T)
#' 
#' alpha <- c(0.03083069, 0.06242601, 0.09120189)
#' lambda <- 0.239385
#' pi <- 0.1453097
#'
#' TT <- 183
#' w <- 28
#' Y <- rzinpoisinar3(TT + w, alpha, lambda, pi, ceiling(TT / 2) + w, delta = 1, burnin = burnin)
#' 
rzinpoisinar3 <- function(n, alpha, lambda, pi, h, delta, burnin = 100) {
  
  q <- length(alpha)
  out <- rep(NA, n + burnin + q)
  out[1:q] <- VGAM::rzipois(q, lambda, pi)
  
  k <- 0
  
  pi1 <- pi1f(delta, lambda, pi)
  lambda1 <- lambda1f(delta, lambda, pi)
  
  for (i in (q + 1):(n + burnin + q)) {
    for (j in 1:q) {
      out[i] <- rbinom(1, out[i - j], alpha[j])
    }
    
    if (i >= (q + 1 + burnin)) {
      k <- k + 1
    }
    
    if (k >= h) {
      out[i] <- out[i] + VGAM::rzipois(1, lambda1, pi1)
    } else {
      out[i] <- out[i] + VGAM::rzipois(1, lambda, pi)
    }
    
    
  }
  
  out[(burnin + q + 1):(n + burnin + q)]
  
}



#' simulate realizations using ARMA(p, q) and one sustained shift
#' 
#' @param n is the length
#' @param phi is the alpha
#' @param theta is the mean of poisson mixture
#' @param sigma2 is the mean of poisson mixture
#' @param h is the proportion of zeros
#' @param delta is the start point of shift
#' @param burnin is the length of the burn-in period
#' @param burnin is the length of the burn-in period
#' @export
#' @examples
#' nsim <- 100
#' burnin <- 100
#' T <- 100
#' q <- 5
#' H <- getHMatMT(T, q)
#' Y <- arima.sim(list(ar = 0.5), n = T)
#' 
#' alpha <- c(0.03083069, 0.06242601, 0.09120189)
#' lambda <- 0.239385
#' pi <- 0.1453097
#'
#' TT <- 183
#' w <- 28
#' Y <- rzinpoisinar3(TT + w, alpha, lambda, pi, ceiling(TT / 2) + w, delta = 1, burnin = burnin)
#' 
rarma <- function(object, n, h, delta, xreg = NULL, nsim = 1000, burnin = 50, lowerbound = 0, const = 0.5) {
  
  ##order <- c(0, 0, 0)
  ##
  ##nar <- sum(object$model$phi != 0)
  ##nma <- sum(object$model$theta != 0)
  ##
  ##if (nar > 0) {
  ##  order[1] <- nar
  ##  phi.vec <- object$model$phi[which(object$model$phi != 0)]
  ##} else {
  ##  phi.vec <- NULL
  ##}
  ##
  ##if (nma > 0) {
  ##  order[3] <- nma
  ##  theta.vec <- object$model$theta[which(object$model$theta != 0)]
  ##} else {
  ##  theta.vec <- NULL
  ##}
  ##
  ##ss <- sigma.mat(100, order = order, phi.vec = phi.vec, theta.vec = theta.vec, sigma2 = object$sigma2, 
  ##                burn.in = burnin)
 
  sim <- matrix(NA, nrow = n, ncol = nsim)
  
  for (i in 1:nsim) {
    sim[, i] <- simulate(object, nsim = n, future = FALSE, xreg = xreg)
  }
  
  fi <- object$fitted
  va <- mean((object$x[1:n] - fi[1:n]) ^ 2)
  
  mu <- rep(0, n)
  mu[h:n] <- mu[h:n] + sqrt(va) * delta
  
  #ts <- simulate(object, nsim = n, future = FALSE, xreg = xreg)
  
  tmpsel <- sample(1:nsim, 1)
  
  ts <- sim[, tmpsel] - const + mu
  
  #innov <- rnorm(n, mu, sqrt(object$sigma2))
  #ts <- simulate(object, nsim = n, future = FALSE, innov = innov, xreg = xreg)

  ts[which(ts < lowerbound)] <- lowerbound
  ts
}


#' simulate realizations using ARMA(p, q) and one sustained shift
#' 
#' @param n is the length
#' @param phi is the alpha
#' @param theta is the mean of poisson mixture
#' @param sigma2 is the mean of poisson mixture
#' @param h is the proportion of zeros
#' @param delta is the start point of shift
#' @param burnin is the length of the burn-in period
#' @param burnin is the length of the burn-in period
#' @export
#' @examples
#' nsim <- 100
#' burnin <- 100
#' T <- 100
#' q <- 5
#' H <- getHMatMT(T, q)
#' Y <- arima.sim(list(ar = 0.5), n = T)
#' 
#' alpha <- c(0.03083069, 0.06242601, 0.09120189)
#' lambda <- 0.239385
#' pi <- 0.1453097
#'
#' TT <- 183
#' w <- 28
#' Y <- rzinpoisinar3(TT + w, alpha, lambda, pi, ceiling(TT / 2) + w, delta = 1, burnin = burnin)
#' 
rarmaAlt <- function(object, n, 
                     h1 = NULL, h2 = NULL, h3 = NULL, 
                     delta1 = 0, delta2 = 0, delta3 = 0, 
                     xreg = NULL, burnin = 50, lowerbound = 0, const = 0.5) {
  
  ##order <- c(0, 0, 0)
  ##
  ##nar <- sum(object$model$phi != 0)
  ##nma <- sum(object$model$theta != 0)
  ##
  ##if (nar > 0) {
  ##  order[1] <- nar
  ##  phi.vec <- object$model$phi[which(object$model$phi != 0)]
  ##} else {
  ##  phi.vec <- NULL
  ##}
  ##
  ##if (nma > 0) {
  ##  order[3] <- nma
  ##  theta.vec <- object$model$theta[which(object$model$theta != 0)]
  ##} else {
  ##  theta.vec <- NULL
  ##}
  ##
  ##ss <- sigma.mat(100, order = order, phi.vec = phi.vec, theta.vec = theta.vec, sigma2 = object$sigma2, 
  ##                burn.in = burnin)
 
  #sim <- matrix(NA, nrow = n, ncol = nsim)
  
  #for (i in 1:nsim) {
  #  sim[, i] <- simulate(object, nsim = n, future = FALSE, xreg = xreg)
  #}
  
  fi <- object$fitted
  va <- mean((object$x[1:n] - fi[1:n]) ^ 2)
  
  n1 <- length(h1)
  n2 <- length(h2)
  n3 <- length(h3)
  
  mu <- rep(0, n)
  
  if (n1 > 0) {
    for (i in 1:n1) {
      mu[h1[i]] <- mu[h1[i]] + sqrt(va) * delta1
    }
  }
  
  if (n2 > 0) {
    for (i in 1:n2) {
      mu[h2[i]:n] <- mu[h2[i]:n] + sqrt(va) * delta2
    }
  }
  
  if (n3 > 0) {
    for (i in 1:n3) {
      mu[h3[i]:n] <- mu[h3[i]:n] + sqrt(va) * delta3 * (h3[i]:n - h3[i] + 1)
    }
  }
  
  
  #ts <- simulate(object, nsim = n, future = FALSE, xreg = xreg)
  
  #tmpsel <- sample(1:nsim, 1)
  
  ts <- simulate(object, nsim = n, future = FALSE, xreg = xreg[1:n, ]) - const + mu
  
  #innov <- rnorm(n, mu, sqrt(object$sigma2))
  #ts <- simulate(object, nsim = n, future = FALSE, innov = innov, xreg = xreg)

  ts[which(ts < lowerbound)] <- lowerbound
  ts
}

#' simulate realizations using ARMA(p, q) and one sustained shift
#'
#' @param n is the length
#' @param phi is the alpha
#' @param theta is the mean of poisson mixture
#' @param sigma2 is the mean of poisson mixture
#' @param h is the proportion of zeros
#' @param delta is the start point of shift
#' @param burnin is the length of the burn-in period
#' @param burnin is the length of the burn-in period
#' @export
#' @examples
#' nsim <- 100
#' burnin <- 100
#' T <- 100
#' q <- 5
#' H <- getHMatMT(T, q)
#' Y <- arima.sim(list(ar = 0.5), n = T)
#'
#' alpha <- c(0.03083069, 0.06242601, 0.09120189)
#' lambda <- 0.239385
#' pi <- 0.1453097
#'
#' TT <- 183
#' w <- 28
#' Y <- rzinpoisinar3(TT + w, alpha, lambda, pi, ceiling(TT / 2) + w, delta = 1, burnin = burnin)
#' 
readrsp <- function(x, lmin = 14, alpha = 0.05) {
  res <- dfphase1::rsp(x, lmin = lmin, alpha = alpha)
  sigmean <- res$p[1] < 0.05
  sigvar <- res$p[2] < 0.05
  sig <- (sigmean + sigvar) > 0
  dd <- diff(res$fit)
  
  dd1 <- 0
  dd2 <- 0
  dd3 <- 0
  
  if (length(which(dd[, 1] > 0)) > 0) {
    dd1 <- which(dd[, 1] > 0) + 1
  }
  
  if (length(which(dd[, 2] > 0)) > 0) {
    dd2 <- which(dd[, 2] > 0) + 1
  }
  
  if ((length(which(dd[, 1] > 0)) + length(which(dd[, 2] > 0))) > 0) {
    dd3 <- sort(unique(c(dd1, dd2)))
  }
  
  out <- list("sigmean" = sigmean, "sigvar" = sigvar, "sig" = sig, 
              "outliermean" = dd1, "outliervar" = dd2, "outlier" = dd3)
  out
}


#' simulate realizations using ARMA(p, q) and one sustained shift
#' 
#' @param n is the length
#' @param phi is the alpha
#' @param theta is the mean of poisson mixture
#' @param sigma2 is the mean of poisson mixture
#' @param h is the proportion of zeros
#' @param delta is the start point of shift
#' @param burnin is the length of the burn-in period
#' @param burnin is the length of the burn-in period
#' @export
#' @examples
#' nsim <- 100
#' burnin <- 100
#' T <- 100
#' q <- 5
#' H <- getHMatMT(T, q)
#' Y <- arima.sim(list(ar = 0.5), n = T)
#' 
#' alpha <- c(0.03083069, 0.06242601, 0.09120189)
#' lambda <- 0.239385
#' pi <- 0.1453097
#'
#' TT <- 183
#' w <- 28
#' Y <- rzinpoisinar3(TT + w, alpha, lambda, pi, ceiling(TT / 2) + w, delta = 1, burnin = burnin)
#' 
rar <- function(ar, n, h, delta, burnin = 50, 
                dist.innov = "norm") {
  
  
  
  if (dist.innov == "norm") {
    innov <- rnorm(n + burnin)
    ss <- sigma.mat(100, order = c(1, 0, 0), phi.vec = ar, theta.vec = NULL, sigma2 = 1, 
                  burn.in = burnin)
    mu <- rep(0, n)
    mu[h:n] <- mu[h:n] + sqrt(ss$gamma0) * delta
  } else if (dist.innov == "gamma") {
    innov <- (rgamma(n + burnin, shape = 1, scale = 1))
    ss <- sigma.mat(100, order = c(1, 0, 0), phi.vec = ar, theta.vec = NULL, 
                    sigma2 = 5 * 25, burn.in = burnin)
    mu <- rep(0, n)
    mu[h:n] <- mu[h:n] + sqrt(ss$gamma0) * delta
  } else if (dist.innov == "poisson") {
    innov <- rpois(n + burnin, 0.6)
    ss <- sigma.mat(100, order = c(1, 0, 0), phi.vec = ar, theta.vec = NULL, 
                    sigma2 = 0.6, burn.in = burnin)
    mu <- rep(0, n)
    mu[h:n] <- mu[h:n] + sqrt(ss$gamma0) * delta
  }
  
  ts <- arima.sim(list(order = c(1, 0, 0), ar = ar), n = n, n.start = burnin, 
                  start.innov = innov[1:burnin], innov = innov[(burnin + 1):(burnin + n)])
  
  ts <- ts + mu
  
  if (dist.innov == "poisson") {
    ts <- round(ts)
  }
  
  #innov <- rnorm(n, mu, sqrt(object$sigma2))
  #ts <- simulate(object, nsim = n, future = FALSE, innov = innov, xreg = xreg)

  ts
}


#' simulate realizations using ARMA(p, q) and one sustained shift
#' 
#' @param n is the length
#' @param phi is the alpha
#' @param theta is the mean of poisson mixture
#' @param sigma2 is the mean of poisson mixture
#' @param h is the proportion of zeros
#' @param delta is the start point of shift
#' @param burnin is the length of the burn-in period
#' @param burnin is the length of the burn-in period
#' @export
#' @examples
#' nsim <- 100
#' burnin <- 100
#' T <- 100
#' q <- 5
#' H <- getHMatMT(T, q)
#' Y <- arima.sim(list(ar = 0.5), n = T)
#' 
#' alpha <- c(0.03083069, 0.06242601, 0.09120189)
#' lambda <- 0.239385
#' pi <- 0.1453097
#'
#' TT <- 183
#' w <- 28
#' Y <- rzinpoisinar3(TT + w, alpha, lambda, pi, ceiling(TT / 2) + w, delta = 1, burnin = burnin)
#' 
rblasso <- function(n, t, delta, model, leftcensoring = 1, rounding = 0, nsim = 1000) {
  
  m <- length(model$Y)
  p <- dim(model$Phi)[1]
  ss <- dim(model$Mu)[2]
  sim <- matrix(NA, nrow = m, ncol = nsim)
  sim[1:p, ] <- model$Y[1:p]
    
  for (i in 1:nsim) {
    tmpsel <- sample(1:ss, 1)
    sim[(p + 1):m, i] <- BayesianLASSOMonitoring::simYph2(0, model$Y, model$Z[, tmpsel], model$Phi[, tmpsel],
                                 model$Mu[, tmpsel], model$sigma2[tmpsel], 1, model$theta[tmpsel], 
                                 leftcensoring, rounding, 1e-6, 1)

  }
  
  rm <- apply(sim, 1, median)
  vv <- mean((model$Y[-c(1:p)] - rm[-c(1:p)]) ^ 2)
  
  tmpsel <- sample(1:nsim, 1)
  out <- sim[1:n, tmpsel]
  out[t:n] <- out[t:n] + sqrt(vv) * delta
  
  if (leftcensoring == 1) {
    out[out < 0] <- 0
  }
  
  if (rounding == 1) {
    out <- round(out)
  }
  
  out
  
} 


#' simulate realizations using ARMA(p, q) and one sustained shift
#' 
#' @param n is the length
#' @param phi is the alpha
#' @param theta is the mean of poisson mixture
#' @param sigma2 is the mean of poisson mixture
#' @param h is the proportion of zeros
#' @param delta is the start point of shift
#' @param burnin is the length of the burn-in period
#' @param burnin is the length of the burn-in period
#' @export
#' @examples
#' nsim <- 100
#' burnin <- 100
#' T <- 100
#' q <- 5
#' H <- getHMatMT(T, q)
#' Y <- arima.sim(list(ar = 0.5), n = T)
#' 
#' alpha <- c(0.03083069, 0.06242601, 0.09120189)
#' lambda <- 0.239385
#' pi <- 0.1453097
#'
#' TT <- 183
#' w <- 28
#' Y <- rzinpoisinar3(TT + w, alpha, lambda, pi, ceiling(TT / 2) + w, delta = 1, burnin = burnin)
#' 
rblassoAlt <- function(model, n, 
                    h1 = NULL, h2 = NULL, h3 = NULL, 
                    delta1 = 0, delta2 = 0, delta3 = 0, 
                    leftcensoring = 1, rounding = 0, nsim = 1000) {
  
  m <- length(model$Y)
  p <- dim(model$Phi)[1]
  ss <- dim(model$Mu)[2]
  sim <- matrix(NA, nrow = m, ncol = nsim)
  sim[1:p, ] <- model$Y[1:p]
    
  for (i in 1:nsim) {
    tmpsel <- sample(1:ss, 1)
    sim[(p + 1):m, i] <- BayesianLASSOMonitoring::simYph2(0, model$Y, model$Z[, tmpsel], model$Phi[, tmpsel],
                                 model$Mu[, tmpsel], model$sigma2[tmpsel], 1, model$theta[tmpsel], 
                                 leftcensoring, rounding, 1e-6, 1)

  }
  
  rm <- apply(sim, 1, median)
  vv <- mean((model$Y[-c(1:p)] - rm[-c(1:p)]) ^ 2)
  
  tmpsel <- sample(1:nsim, 1)
  out <- sim[1:n, tmpsel]
  
  n1 <- length(h1)
  n2 <- length(h2)
  n3 <- length(h3)
  
  mu <- rep(0, n)
  
  if (n1 > 0) {
    for (i in 1:n1) {
      mu[h1[i]] <- mu[h1[i]] + sqrt(vv) * delta1
    }
  }
  
  if (n2 > 0) {
    for (i in 1:n2) {
      mu[h2[i]:n] <- mu[h2[i]:n] + sqrt(vv) * delta2
    }
  }
  
  if (n3 > 0) {
    for (i in 1:n3) {
      mu[h3[i]:n] <- mu[h3[i]:n] + sqrt(vv) * delta3 * (h3[i]:n - h3[i] + 1)
    }
  }
  
  out <- out + mu
  
  if (leftcensoring == 1) {
    out[out < 0] <- 0
  }
  
  if (rounding == 1) {
    out <- round(out)
  }
  
  out
  
} 

#' simulate realizations using ARMA(p, q) and one sustained shift
#' 
#' @param n is the length
#' @param phi is the alpha
#' @param theta is the mean of poisson mixture
#' @param sigma2 is the mean of poisson mixture
#' @param h is the proportion of zeros
#' @param delta is the start point of shift
#' @param burnin is the length of the burn-in period
#' @param burnin is the length of the burn-in period
#' @export
#' @examples
#' nsim <- 100
#' burnin <- 100
#' T <- 100
#' q <- 5
#' H <- getHMatMT(T, q)
#' Y <- arima.sim(list(ar = 0.5), n = T)
#' 
#' alpha <- c(0.03083069, 0.06242601, 0.09120189)
#' lambda <- 0.239385
#' pi <- 0.1453097
#'
#' TT <- 183
#' w <- 28
#' Y <- rzinpoisinar3(TT + w, alpha, lambda, pi, ceiling(TT / 2) + w, delta = 1, burnin = burnin)
#' 
rblassoAlt1 <- function(model, n, 
                    h1 = NULL, h2 = NULL, h3 = NULL, 
                    delta1 = 0, delta2 = 0, delta3 = 0, 
                    leftcensoring = 1, rounding = 0, nsim = 1000) {
  
  m <- length(model$Y)
  p <- dim(model$Phi)[1]
  ss <- dim(model$Mu)[2]
  sim <- matrix(NA, nrow = m, ncol = nsim)
  sim[1:p, ] <- model$Y[1:p]
    
  for (i in 1:nsim) {
    tmpsel <- sample(1:ss, 1)
    sim[(p + 1):m, i] <- BayesianLASSOMonitoring::simYph2Alt(0, model$Y, model$Z[, tmpsel], model$Phi[, tmpsel],
                                 model$Mu[, tmpsel], model$sigma2[tmpsel], 1, model$theta[tmpsel], 
                                 leftcensoring, rounding, 1e-6, 1)

  }
  
  rm <- apply(sim, 1, median)
  vv <- mean((model$Y[-c(1:p)] - rm[-c(1:p)]) ^ 2)
  
  tmpsel <- sample(1:nsim, 1)
  out <- sim[1:n, tmpsel]
  
  n1 <- length(h1)
  n2 <- length(h2)
  n3 <- length(h3)
  
  mu <- rep(0, n)
  
  if (n1 > 0) {
    for (i in 1:n1) {
      mu[h1[i]] <- mu[h1[i]] + sqrt(vv) * delta1
    }
  }
  
  if (n2 > 0) {
    for (i in 1:n2) {
      mu[h2[i]:n] <- mu[h2[i]:n] + sqrt(vv) * delta2
    }
  }
  
  if (n3 > 0) {
    for (i in 1:n3) {
      mu[h3[i]:n] <- mu[h3[i]:n] + sqrt(vv) * delta3 * (h3[i]:n - h3[i] + 1)
    }
  }
  
  out <- out + mu
  
  if (leftcensoring == 1) {
    out[out < 0] <- 0
  }
  
  if (rounding == 1) {
    out <- round(out)
  }
  
  out
  
} 

#' Caculate the moving averages
#' 
#' gets the moving averages
#' @param Y is the input
#' @param w is the length of moving window
#' 
#' @export
#' @examples
#' alpha <- c(0.03083069, 0.06242601, 0.09120189)
#' lambda <- 0.239385
#' pi <- 0.1453097
#'
#' TT <- 183
#' w <- 28
#' Y <- rzinpoisinar3(TT + w, alpha, lambda, pi, ceiling(TT / 2) + w, delta = 1, burnin = 100)
#' ma <- movaver(Y, w)
movaver <- function(Y, w = 5){filter(Y, rep(1 / w, w), sides = 1)}

#' Transform the data
#' 
#' gets the transformed input
#' @param Y is the input
#' @param log is the flag triggering the log transformation
#' @param const is the constant added to the input during the log transformation
#' @param sta is the flag triggering the standardization
#' 
#' @export
#' @examples
#' alpha <- c(0.03083069, 0.06242601, 0.09120189)
#' lambda <- 0.239385
#' pi <- 0.1453097
#'
#' TT <- 183
#' w <- 28
#' Y <- rzinpoisinar3(TT + w, alpha, lambda, pi, ceiling(TT / 2) + w, delta = 1, burnin = 100)
#' ma <- movaver(Y, w)
#' ma.tr <- trans(ma, TRUE, 1, TRUE)
trans <- function(Y, log = TRUE, const = 1, sta = TRUE, meanY = NULL, sdY = NULL){
  out <- Y
  if (log == TRUE) {
    out <- log(out + const)
  }
  if (sta == TRUE) {
    meanY <- ifelse(!is.null(meanY), meanY, mean(out))
    sdY <- ifelse(!is.null(sdY), sdY, sd(out))
    out <- (out - meanY) / sdY
  }
  
  list(
    "Y" = out,
    "meanY" = meanY,
    "sdY" = sdY
  )
}

#' Back-transform the data
#' 
#' gets the back-transformed input
#' @param Y is the input
#' @param log is the flag triggering the log transformation
#' @param const is the constant added to the input during the log transformation
#' @param sta is the flag triggering the standardization
#' @param meanY is the mean of the original Y
#' @param sdY is the standard deviation of the original Y
#' 
#' @export
#' @examples
#' alpha <- c(0.03083069, 0.06242601, 0.09120189)
#' lambda <- 0.239385
#' pi <- 0.1453097
#'
#' TT <- 183
#' w <- 28
#' Y <- rzinpoisinar3(TT + w, alpha, lambda, pi, ceiling(TT / 2) + w, delta = 1, burnin = 100)
#' ma <- movaver(Y, w)
#' ma.tr <- trans(ma, TRUE, 1, TRUE)
#' ma.batr <- backtrans(ma.tr$Y, TRUE, 1, TRUE, ma.tr$meanY, ma.tr$sdY)
backtrans <- function(Y, log = TRUE, const = 1, sta = TRUE, meanY = 0, sdY = 1){
  out <- Y
  
  if (sta == TRUE) {
    out <- out * sdY + meanY
  }
  
  if (log == TRUE) {
    out <- exp(out) - const
  }
  
  out
}

#' simulate realizations using INAR(3) with zero-inflated Poisson innovation and one sustained shift
#' 
#' @param Ph1BayesianLASSO.model is the length
#' @param log is the log
#' @param const is the constant
#' @param sta is the sta
#' @export
#' @examples
#' nsim <- 100
#' burnin <- 100
#' T <- 100
#' q <- 5
#' H <- getHMatMT(T, q)
#' Y <- arima.sim(list(ar = 0.5), n = T)
#' 
#' alpha <- c(0.03083069, 0.06242601, 0.09120189)
#' lambda <- 0.239385
#' pi <- 0.1453097
#'
#' TT <- 183
#' w <- 28
#' Y <- rzinpoisinar3(TT + w, alpha, lambda, pi, ceiling(TT / 2) + w, delta = 1, burnin = burnin)
#' 
RMSE.ph1 <- function(Ph1BayesianLASSO.model, log = TRUE, const = 1, sta = TRUE, lowerbound = 0) {
  nsim <- dim(Ph1BayesianLASSO.model$Phi)[2]
  q <- dim(Ph1BayesianLASSO.model$Phi)[1]
  TT <- length(Ph1BayesianLASSO.model$Y.tr)
  
  tmpfit.tr <- rep(NA, TT)
  tmpfit.ma <- rep(NA, TT)
  
  tmpresi.tr <- rep(NA, TT)
  tmpresi.ma <- rep(NA, TT)
  
  RMSE.tr <- rep(NA, nsim)
  RMSE.ma <- rep(NA, nsim)
  
  fit.tr <- matrix(NA, nrow = TT, ncol = nsim)
  fit.ma <- matrix(NA, nrow = TT, ncol = nsim)
  
  resi.tr <- matrix(NA, nrow = TT, ncol = nsim)
  resi.ma <- matrix(NA, nrow = TT, ncol = nsim)
  
  Y <- Ph1BayesianLASSO.model$Y
  Y.tr <- Ph1BayesianLASSO.model$Y.tr
  Y.ma <- Ph1BayesianLASSO.model$Y.ma
  meanY <- Ph1BayesianLASSO.model$meanY
  sdY <- Ph1BayesianLASSO.model$sdY
  
  X <- Ph1BayesianLASSO.model$X
  H <- Ph1BayesianLASSO.model$H
  
  for (j in 1:nsim) {
    tmpmuq <- Ph1BayesianLASSO.model$muq[j]
    tmpPhi <- Ph1BayesianLASSO.model$Phi[, j]
    tmpV <- Y.tr - tmpmuq
    if (!is.null(X)) {
      tmpBeta <- Ph1BayesianLASSO.model$Beta[, j]
      tmpKappa <- Ph1BayesianLASSO.model$Kappa[, j]
      tmpV <- tmpV - X %*% (tmpBeta * tmpKappa) 
    }
    if (!is.null(H)) {
      tmpGamma <- Ph1BayesianLASSO.model$Gamma[, j]
      tmpTau <- Ph1BayesianLASSO.model$Tau[, j]
      tmpV <- tmpV - H %*% (tmpGamma * tmpTau)
    }
    for (i in (q + 1):TT) {
      tmpfit.tr[i] <- tmpmuq + tmpV[(i - 1):(i - q)] %*% tmpPhi
      
      if (!is.null(X)) {
        tmpfit.tr[i] <- tmpfit.tr[i] + X[i, ] %*% (tmpBeta * tmpKappa)
      }
      
      
      if (!is.null(H)) {
        tmpfit.tr[i] <- tmpfit.tr[i] + H[i, ] %*% (tmpGamma * tmpTau)
      }
      
      tmpresi.tr[i] <- Y.tr[i] - tmpfit.tr[i]
      tmpfit.ma[i] <- backtrans(tmpfit.tr[i], log, const, sta, meanY, sdY)
      
      tmpfit.ma[i] <- ifelse(tmpfit.ma[i] < lowerbound, lowerbound, tmpfit.ma[i])
      tmpresi.ma[i] <- Y.ma[i] - tmpfit.ma[i]
      
    }
    
    fit.tr[, j] <- tmpfit.tr
    fit.ma[, j] <- tmpfit.ma
    
    resi.tr[, j] <- tmpresi.tr
    resi.ma[, j] <- tmpresi.ma

    tmpresi.tr <- tmpresi.tr[(q + 1):TT]
    RMSE.tr[j] <- sqrt(t(tmpresi.tr) %*% tmpresi.tr / (TT - q))
    
    tmpresi.ma <- tmpresi.ma[(q + 1):TT]
    RMSE.ma[j] <- sqrt(t(tmpresi.ma) %*% tmpresi.ma / (TT - q))
    
  }
  
  list("RMSE.tr" = RMSE.tr, "RMSE.ma" = RMSE.ma, 
       "fit.tr" = fit.tr[-c(1:q), ], "fit.ma" = fit.ma[-c(1:q), ], 
       "resi.tr" = resi.tr[-c(1:q), ], "resi.ma" = resi.ma[-c(1:q), ])

}

#' simulate realizations using INAR(3) with zero-inflated Poisson innovation and one sustained shift
#' 
#' @param Ph1BayesianLASSO.model is the length
#' @param log is the log
#' @param const is the constant
#' @param sta is the sta
#' @export
#' @examples
#' nsim <- 100
#' burnin <- 100
#' T <- 100
#' q <- 5
#' H <- getHMatMT(T, q)
#' Y <- arima.sim(list(ar = 0.5), n = T)
#' 
#' alpha <- c(0.03083069, 0.06242601, 0.09120189)
#' lambda <- 0.239385
#' pi <- 0.1453097
#'
#' TT <- 183
#' w <- 28
#' Y <- rzinpoisinar3(TT + w, alpha, lambda, pi, ceiling(TT / 2) + w, delta = 1, burnin = burnin)
#' 
RMSE.ph2 <- function(Y, Ph1BayesianLASSO.model, X = NULL, H = NULL, 
                     log = TRUE, const = 1, sta = TRUE, lowerbound = 0) {
  
  TT <- length(Y)
  YY <- c(Ph1BayesianLASSO.model$Y, Y)
  nn <- length(YY)
  q <- dim(Ph1BayesianLASSO.model$Phi)[1]
  
  Y.ma <- movaver(YY, Ph1BayesianLASSO.model$w)[(nn - TT - q + 1):nn]
  tt <- length(Y.ma)
  
  meanY <- Ph1BayesianLASSO.model$meanY
  sdY <- Ph1BayesianLASSO.model$sdY
  
  Y.tr <- trans(Y.ma, log, const, sta, meanY, sdY)$Y
  
  nsim <- dim(Ph1BayesianLASSO.model$Phi)[2]

  
  tmpfit.tr <- rep(NA, tt)
  tmpfit.ma <- rep(NA, tt)
  tmpresi.tr <- rep(NA, tt)
  tmpresi.ma <- rep(NA, tt)
  RMSE.tr <- rep(NA, nsim)
  RMSE.ma <- rep(NA, nsim)
  
  if (!is.null(X)) {
    X <- as.matrix(rbind(Ph1BayesianLASSO.model$X, X)[(nn - TT - q + 1):nn, ])
  }
  if (!is.null(H)) {
    H <- as.matrix(rbind(Ph1BayesianLASSO.model$H, H)[(nn - TT - q + 1):nn, ])
  }
  
  for (j in 1:nsim) {
    tmpmuq <- Ph1BayesianLASSO.model$muq[j]
    tmpPhi <- Ph1BayesianLASSO.model$Phi[, j]
    tmpV <- Y.tr - tmpmuq
    if (!is.null(X)) {
      tmpBeta <- Ph1BayesianLASSO.model$Beta[, j]
      tmpKappa <- Ph1BayesianLASSO.model$Kappa[, j]
      tmpV <- tmpV - X %*% (tmpBeta * tmpKappa) 
    }
    if (!is.null(H)) {
      tmpGamma <- Ph1BayesianLASSO.model$Gamma[, j]
      tmpTau <- Ph1BayesianLASSO.model$Tau[, j]
      tmpV <- tmpV - H %*% (tmpGamma * tmpTau)
    }
    
    for (i in (q + 1):tt) {
      tmpfit.tr[i] <- tmpmuq + tmpV[(i - 1):(i - q)] %*% tmpPhi
      
      if (!is.null(X)) {
        tmpfit.tr[i] <- tmpfit.tr[i] + X[i, ] %*% (tmpBeta * tmpKappa)
      }
      if (!is.null(H)) {
        tmpfit.tr[i] <- tmpfit.tr[i] + H[i, ] %*% (tmpGamma * tmpTau)
      }
      
      tmpresi.tr[i] <- Y.tr[i] - tmpfit.tr[i]
      tmpfit.ma[i] <- backtrans(tmpfit.tr[i], log, const, sta, meanY, sdY)
      tmpfit.ma[i] <- ifelse(tmpfit.ma[i] < lowerbound, lowerbound, tmpfit.ma[i])
      tmpresi.ma[i] <- Y.ma[i] - tmpfit.ma[i]
    }
    tmpresi.tr <- tmpresi.tr[(q + 1):tt]
    RMSE.tr[j] <- sqrt(t(tmpresi.tr) %*% tmpresi.tr / (tt - q))
    
    tmpresi.ma <- tmpresi.ma[(q + 1):tt]
    RMSE.ma[j] <- sqrt(t(tmpresi.ma) %*% tmpresi.ma / (tt - q))
  }
  
  list("RMSE.tr" = RMSE.tr, "RMSE.ma" = RMSE.ma)
  
}
