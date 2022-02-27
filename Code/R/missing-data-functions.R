library(fastfactoranalysis)

genMask <- function(Y, nmiss, type="MCAR"){
  n <- nrow(Y)
  p <- ncol(Y)
  if(type == "MCAR"){
    s <- sample(n, nmiss)
    M <- matrix(0, n, p)
    M[s, 1] <- 1
  } else if(type == "trunc"){
    M <- matrix(0,n, p)
    M[Y[,1] > quantile(Y[,1], (1-nmiss/n)),1] <- 1
  } else {
    stop("imputation type unknown")
  }
  M>0
}

impute <- function(Y, M){
  Ymeans <- colMeans(Y)
  Yimp <- sapply(1:ncol(Y), function(j){
    Yj <- Y[,j]
    Yj[M[,j]] <- Ymeans[j]
    Yj
  })
  Yimp
}


genY <- function(n, p, type="IID", ...){
  if(type=="IID"){
    args <- list(...)
    if(is.null(args$mu)) mu <- rep(0, p) else mu <- args$mu
    if(is.null(args$psi)) psi <- rep(1, p) else psi <- args$psi
    Y <- sapply(1:p, function(j){
      rnorm(n, mu[j], psi[j])
    })
    args <- list(mu=mu)
  }
  if(type=="FA"){
    args <- list(...)
    if(is.null(args$q))stop("Please specify q.")
    if(is.null(args$psi)) psi <- rep(1, p) else psi <- args$psi
    if(is.null(args$mu)) mu <- rep(0, p) else mu <- args$mu
    if(is.null(args$A)) A <- matrix(rnorm(p*q), p, q) else A <- args$A
    if(is.null(args$Z)) Z <- matrix(rnorm(n*q), n, q) else Z <- args$Z
    linpar <- t(A %*% t(Z) + mu)
    Y <- sapply(1:p, function(j){
      rnorm(n, linpar[,j], psi[j])
    })
    args <- list(mu=mu, A=A, linpar=linpar, psi=psi, Z=Z)
  }
  list(Y=Y, args=args)
}

theta_to_par <- function(theta, p, q){
  A <- matrix(theta[1:(p*q)], p, q)
  psi <- theta[(p*q+1):(p*q+p)]
  mu  <- theta[(p*q+p + 1): length(theta)]
  list(A=A, psi=psi, mu=mu)
}

par_to_theta <- function(A, psi, mu){
  c(as.vector(A), psi, mu)
}

# dat is a list with data elements that will be the first argument of estimator
# estimator must take a list named dat as first parameter, the rest in est.args.
ib <- function(dat, estimator, generator, theta.init=NULL, H=1, step=0.1, maxit=100, gen.args=list(), est.args=list(), seed=12451, verbose=F, save=T,
               momentum=0, tol=1e-6, method="SA", method.args=list()){
  target <- theta <- do.call(estimator, c(list(dat=dat), est.args))
  if(!is.null(theta.init)) theta <- theta.init

  if(save){
    hist <- matrix(0, maxit, length(target))
    hist[1,] <- theta
  }
  theta.avg.old <- theta

  mom <- rep(0, length(theta))
  delta <- rep(0, length(theta))
  learning.rate <- 1

  if(method=="SA"){
    if(is.null(method.args$a)) method.args$a <- 40
    if(is.null(method.args$A)) method.args$A <- .05*maxit
    if(is.null(method.args$alpha)) method.args$alpha <- .8
    if(is.null(method.args$b)) method.args$b <- 1

  }


  for(i in 2:maxit){
    delta.old <- delta
    theta.old <- theta

    exp <- lapply(1:H, function(h){
      if(method=="SP") set.seed(seed+h)
      dat.gen <- do.call(generator, c(list(theta=theta), gen.args))
      est <- do.call(estimator, c(list(dat=dat.gen), est.args))
      est
    })
    exp <- Reduce("+", exp)/length(exp)

    # update
    if(method=="SA") learning.rate <- with(method.args, a/(b*i + 1 + A)**alpha)
    cat("\n", i, " : ", learning.rate)
    # compute the step
    delta <- learning.rate * step* (target-exp)
    theta <- theta + delta + mom
    mom <- (mom + delta)*momentum

    theta[is.na(theta)] <- -Inf
    if(save) hist[i,] <- theta

    if(i %% 10 == 0){
      if(save){
        theta.avg <- avg.last.iter(hist[1:i,], i/3)
        crit <- sqrt(mean((theta.avg.old-theta.avg)**2))
        theta.avg.old <- theta.avg
      } else {
        crit <- sqrt(mean((theta.old-theta)**2))
      }
      cat(" - ", crit)
      if(crit < tol){
        if(save) hist <- hist[1:i,]
        break()
      }
    }
  }
  # returns the average over the last 10% iterations
  list(theta=theta, target=target, save= if(save) hist else NULL, theta.avg=if(save) avg.last.iter(hist, nrow(hist)/3))
}

# returns the average over  the last iterations

avg.last.iter <- function(theta.hist, last.iter){
  stopifnot(is.matrix(theta.hist))
  last.iter <- min(nrow(theta.hist), last.iter)
  colMeans(theta.hist[(nrow(theta.hist) + 1 - last.iter):nrow(theta.hist),, drop=F])
}

# #########################
# Part 3: ADDING THE MASK AND LEARNING EVERYTHING
# LOW RANK APPROXIMATION WITH FACTOR ANALYSIS

# given the data, get various estimators of the mean...
generator_M_known <- function(theta, n, p, q, nmiss, missType){
  par <- theta_to_par(theta, p, q)
  Z <- matrix(rnorm(n*q), n ,q)
  linpar <- t(par$A %*% t(Z) + par$mu)
  Y <- sapply(1:p, function(j){
    rnorm(n, linpar[,j], par$psi[j])
  })
  M <- genMask(Y, nmiss, missType)
  list(Y=Y, M=M)
}

# Y is a list with Y$Y the full data, Y$M the mask
estimator_M_known <- function(Y, q){
  p <- ncol(Y$Y)
  Yimp <- impute(Y$Y, Y$M)
  fit <- ffa(Yimp, nfactors=q)

  par_to_theta(A=fit$loadings, psi=fit$communalites, mu=fit$means)
}

# pmiss: number of variables with missing values
generator_M_unknown <- function(theta, n, p, q){
  par <- theta_to_par(theta, p, q)
  Z <- matrix(rnorm(n*q), n ,q)
  linpar <- t(par$A %*% t(Z) + par$mu)
  Y <- sapply(1:(p/2), function(j){
    rnorm(n, linpar[,j], par$psi[j])
  })
  M <- sapply((p/2+1):p, function(j){
    rbinom(n, 1, 1/(1+exp(-linpar[,j])))
  })>0
  list(Y=Y, M=M)
}

# Y is a list with Y$Y the full data, Y$M the mask
estimator_M_unknown <- function(Y, q){
  p <- ncol(Y$Y)*2
  no.missing <- c(rep(T, p/2), colMeans(Y$M)!=0)
  Yimp <- impute(Y$Y, Y$M)
  fit <- ffa(cbind(Yimp, Y$M*4)[, no.missing], nfactors=q)

  # Add rows of 0 for the missing variables
  A <- matrix(0, p, q)
  A[no.missing, ] <- fit$loadings
  psi <- rep(.1, p)
  psi[no.missing] <- fit$communalites
  psi[(p/2+1):p] <- .1
  mu <- rep(-Inf, p)
  mu[no.missing] <- fit$means
  par_to_theta(A=A, psi=psi, mu=mu)
}
