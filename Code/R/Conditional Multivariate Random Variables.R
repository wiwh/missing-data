library(mvtnorm)
source("R/missing-data-functions.R")
# Generating Conditional Multivariate Normal Random Variables
# We reproduce results from Josse


# generate random variables from the multivariate normal conditional distribution
#'@param n number of random variables to draw
#'@param X observed values conditioned upon, of size p. X[!condvar] is ignored.
#'@param condvar a vector of size p true / false denoting the conditioning variables
#'@param mu a vector of size p containing the mean of the random variables
#'@param sigma a matrix containing
#'@value a conditioning variable... whateve.
rmnorm.cond <- function(n, x, condvar, mu=rep(0, length(condvar)), sigma = diag(length(condvar))){
  sigma11 = sigma[!condvar, !condvar, drop=F]
  sigma12 = sigma[!condvar, condvar, drop=F]
  sigma22 = sigma[condvar, condvar, drop=F]
  sigma21 = sigma[condvar, !condvar, drop=F]
  # if sigma.inverse is given, we use the schur complement.
  sigma22.inverse <- solve(sigma22)
  mu.cond <- mu[!condvar] + as.vector(sigma12 %*% sigma22.inverse %*% (x[condvar] - mu[condvar]))
  sigma.cond <- sigma11 - sigma12 %*% sigma22.inverse %*% sigma21

  list(x=mvtnorm::rmvnorm(n, mean=mu.cond, sigma=sigma.cond), mu=mu.cond, sigma=sigma.cond)
}

mean.obs <- function(X, M){
  X[M] <- NA
  colMeans(X, na.rm=T)
}

theta_to_args <- function(theta, p){
  mu <- theta[1:p]
  Sigma <- matrix(theta[-(1:p)],p ,p )
  list(mu=mu, Sigma=Sigma)
}

args_to_theta <- function(args){
  c(args$mu, as.vector(args$Sigma))
}

# TODO: write this
#'@param dat: list of (X, M)
estimator <- function(dat, type="meanImpute"){
  if(type=="MAP"){
    args <- theta_to_args(dat$theta, ncol(dat$X))
    Ximp <- t(sapply(1:nrow(dat$X), function(i){
      mu.cond <- args$mu
      if(sum(!dat$M[i,])==0) stop(paste("All variables missing for observation", i))
      if(sum(!dat$M[i,])==length(!dat$M[i,])){
        warning(paste("Conditioning on all variables for observation", i))
      } else {
        mu.cond[dat$M[i,]] <- rmnorm.cond(1,dat$X[i,], condvar=!dat$M[i,], mu=args$mu, sigma=args$Sigma)$mu
      }
      mu.cond
    }))
  }
  if(type=="meanImpute"){
    Ximp <- impute(dat$X, dat$M)
  }
  args <- list(
    mu = colMeans(Ximp),
    Sigma = cov(Ximp)
  )
  args_to_theta(args)
}

generator <- function(theta, X, M, p){
  args <- theta_to_args(theta, p)
  Ximp <- t(sapply(1:nrow(X), function(i){
    x <- X[i,]
    if(sum(!M[i,])==0) stop(paste("All variables missing for observation", i))
    if(sum(!M[i,])==length(!M[i,])){
      warning(paste("Conditioning on all variables for observation", i))
    } else {
      x[M[i,]] <- rmnorm.cond(1,X[i,], condvar=!M[i,], mu=args$mu, sigma=args$Sigma)$x
    }
    x
  }))
  list(X=Ximp, M=M, theta=theta)
}

generator.unconditional <- function(theta, X, M, p){
  args <- theta_to_args(theta, p)
  X <- rmvnorm(n, mean=args$mu, sigma = args$Sigma)
  list(X=X, M=M, theta=theta)
}

# simus
n <- 1000
p <- 5
C <- matrix(c(
  1, 0.8, 0, 0, 0,
  0.8, 1, 0, 0, 0,
  0, 0, 1, 0.3, 0.6,
  0, 0, 0.3, 1, 0.7,
  0, 0, 0.6, 0.7, 1
), ncol=p, byrow=T)

sigma <- 1:p
mu <- 1:p

Sigma <- diag(sigma) %*% C %*% diag(sigma)
args0 <- list(mu=mu, Sigma=Sigma)

set.seed(213411)
X <- rmvnorm(n, mean=args0$mu, sigma=args0$Sigma)
M <- matrix(runif(n*p)<.1, n, p)
dat <- list(X=X, M=M, theta=NULL)
theta0 <- args_to_theta(args0)


X2 <- generator(theta0, X, M=M, p)$X
dat2 <- list(X=X2, M=M, theta=NULL)

estimator(dat2, type="meanImpute")
estimator(dat, type="meanImpute")


gen.args <- list(X=X, M=M, p=p)
est.args <- list(type="meanImpute")


ib.theta <- ib(dat, estimator, generator.unconditional,
            gen.args = gen.args,
            est.args = est.args,
            save = T, maxit = 5, H=10, step=1)

# ib.theta <- ib(dat, estimator, generator,
#             gen.args = gen.args,
#             est.args = est.args,
#             save = T, maxit = 100, H=10)

ib.args <- theta_to_args(ib.theta$theta, p)

ts.plot(ib.theta$save)


n <- 10000
X <- rmvnorm(n, mean=mu, sigma=Sigma)
M <- matrix(runif(n*p),n, p)<.5

X1 <- t(sapply(1:n, function(i){
  x <- X[i,]
  Msum <- sum(M[i,])
  if(!(Msum %in% c(0, p))){
    x[!M[i,]] <- rmnorm.cond(1, X[i,], M[i,], mu=mu, sigma=Sigma)$x
  }
  x
  })
)
X2 <- t(sapply(1:n, function(i){
  x <- X[i,]
  if(!(sum(!M[i,]) %in% c(0, p))){
    x[M[i,]] <- rmnorm.cond(1, X[i,], !M[i,], mu=mu, sigma=Sigma)$x
  }
  x
  })
)

plot(cov(X), cov(X1)); abline(0,1,col=2)
points(cov(X), cov(X2),col=2)
