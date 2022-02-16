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

# Y is a list with Y$Y the complete data, Y$M the mask.
# estimator must take a list named Y as first parameter, the rest in est.args.
# generator must take theta, a vector, as first parameter, the rest in gen.args
ib <- function(Y, estimator, generator, theta.init=NULL, H=1, step=0.1, maxit=100, gen.args=list(), est.args=list(), seed=12451, verbose=F, save=F){
  target <- theta <- do.call(estimator, c(list(Y=Y), est.args))
  if(!is.null(theta.init)) theta <- theta.init
  if(save) hist <- matrix(0, maxit, length(target))
  for(i in 1:maxit){
    exp <- lapply(1:H, function(h){
      set.seed(seed+h)
      Y <- do.call(generator, c(list(theta=theta), gen.args))
      est <- do.call(estimator, c(list(Y=Y), est.args))
      est
    })
    exp <- Reduce("+", exp)/length(exp)
    theta <- theta + step *(target - exp)
    theta[is.na(theta)] <- -Inf
    cat(i, "\n", theta[c(1,6)])
    if(save) hist[i,] <- theta
  }
  list(theta=theta, target=target, save= if(save) hist else NULL)
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
