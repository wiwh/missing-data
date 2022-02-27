
source("R/missing-data-functions.R")

args_to_theta <- function(args){
  c(args$beta, args$mu, as.vector(args$Sigma))
}

theta_to_args <- function(theta, p){
  beta <- theta[1:(p+1)]
  mu <- theta[(p+2):(2*p+1)]
  Sigma <- matrix(theta[-(1:(2*p+1))], p, p)
  list(beta=beta, mu=mu, Sigma=Sigma)
}


generator <- function(theta, n, p){
  args <- theta_to_args(theta, p)
  # Generate X
  X <- mvtnorm::rmvnorm(n, mean=args$mu, sigma=args$Sigma)
  linpar <- cbind(1, X) %*% args$beta
  mu <- 1/(1+exp(-linpar))
  Y <- rbinom(nrow(X),1, mu)
  list(Y=Y, X=X, linpar=linpar, mu=mu)
}

estimator.glm <- function(dat, M){
  dat$X <- impute(dat$X, M)
  args <- list(
    beta = glm(Y~ 1 + X, data=dat, family="binomial")$coef, # intercept gets included here
    mu = colMeans(dat$X),
    Sigma = cov(dat$X)
  )
  args_to_theta(args)
}

estimator.glmnet <- function(dat, M, lambda=.1){
  dat$X <- impute(dat$X, M)
  fit.glmnet = glmnet(dat$X, dat$Y, family="binomial", alpha=0, lambda=lambda) # intercept gets included here
  beta <- c(fit.glmnet$a0, as.vector(fit.glmnet$beta))
  args <- list(
    beta=beta,
    mu = colMeans(dat$X),
    Sigma = cov(dat$X)
  )
  args_to_theta(args)
}

estimator.lm <- function(dat, M){
  dat$X <- impute(dat$X, M)
  args <- list(
    beta = lm(Y~ 1 + X, data=dat)$coef, # intercept gets included here
    mu = colMeans(dat$X),
    Sigma = cov(dat$X)
  )
  args_to_theta(args)
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

mu <- 1:p

sigma <- runif(5)
Sigma <- diag(sigma) %*% C %*% diag(sigma)
beta <- c(-.2, .5, -.3, 1, 0, -.6)  # the first element is the intercept
args0 <- list(beta=beta, mu=mu, Sigma=Sigma)
theta0 <- args_to_theta(args0)
prop.miss <- .1
M <- matrix(runif(n*p)<prop.miss, n, p)


simpar <- list(
  n=n,
  p=p,
  Sigma=Sigma,
  beta=beta,
  args0=args0,
  theta0=theta0,
  M=M
)

set.seed(52500)

dat0 <- generator(theta0, n, p)
glm.est <- glm(Y~ X,data=dat0, family="binomial")$coef

plot(estimator.glm(dat0, M), theta0)

stepsize <- rep(1, length(theta0))
stepsize[1:(p+1)] <- 1


fit.ib <- ib(dat0, estimator.lm, generator, method="SP", method.args=method.args, H=100,  est.args=list(M=M), gen.args = list(n=n, p=p), save=T, maxit=50, step=stepsize, momentum=0.6)
par(mfrow=c(2,1))
ts.plot(fit.ib$save)
ts.plot(fit.ib$save[,1:6])
par(mfrow=c(1,1))


maxit.sa <- 200
method.args <- list(a=20, A=maxit.sa*.05, alpha=1.2, b=1)
fit.ib.sa <- ib(dat0, estimator.lm, generator,  method="SA", method.args=method.args, H=1,  est.args=list(M=M), gen.args = list(n=n, p=p), save=T, maxit=maxit.sa, step=stepsize, momentum=0.6)

par(mfrow=c(2,1))
ts.plot(fit.ib.sa$save)
ts.plot(fit.ib.sa$save[,1:6])
par(mfrow=c(1,1))

est.glm <- glm(Y~X, data=dat0, family=binomial)$coef
X.impute <- impute(dat0$X, M)
est.glm.impute <- glm(dat0$Y ~X.impute, family=binomial)$coef

args.sa.avg <- theta_to_args(fit.ib$theta.avg, p)
args <- theta_to_args(fit.ib$theta, p)
args.sa <- theta_to_args(fit.ib.sa$theta, p)

plot(args0$beta, args$beta); abline(0,1,col=2)
points(args0$beta, args.sa$beta, col=3, pch=2)
points(args0$beta, args.sa.avg$beta, col=2, pch=2)
points(args0$beta, est.glm, col=3)
points(args0$beta, est.glm.impute, col=3, pch=3)
