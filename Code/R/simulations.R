source("./R/adding-the-mask_functions.R")

# only works when p =3...
get_theta_init <- function(par.init){
  par.init$A <- rbind(par.init$A, matrix(c(5, 0, 0)))
  par.init$psi <- rep(par.init$psi, 2)
  par.init$mu <- c(par.init$mu, c(0, -Inf, -Inf))
  theta.init <- do.call(par_to_theta, par.init)
  theta.init
}

n <- 1000
p <- 3
q <- 1
nmiss <- .5*n
missType <- "MCAR"
missType <- "trunc"
set.seed(10243)
A <- matrix(rnorm(p*q), p, q)*1.5
mu <- c(c(0, -1, 1), rnorm(p-3))
psi <- rep(0.5, p)
A[1:3] <- c(1, 1, -1)


gen_dat <- function(seed){
  set.seed(seed)
  dat0 <- genY(n, p, q=q, type="FA", A=A, psi=psi, mu=mu)
  dat0$M <- genMask(dat0$Y, nmiss = nmiss, type=missType)
  dat0
}

# 1 : Generate an example of missing data
if(0){
dat0 <- gen_dat(132312)
y0 <- dat0$Y
colnames(y0) <- c("Y1", "Y2", "Y3")
y0[dat0$M] <- NA
View(round(y0, 2))

miss <- dat0$M *1
colnames(miss) <- c("M1", "M2", "M3")
View(miss)
}


# 2: experiment!

# We need 3 different estimators:

# naive: impute and take avg
est_naive_impute <- function(dat){
  Yimp <- impute(dat$Y, dat$M)
  colMeans(Yimp)
}

# bias-corrected, known M
est_M_known <- function(dat){
  theta.init <- par_to_theta(dat$args$A, dat$args$psi, dat$args$mu)
  fit  <- ib(dat[c("Y", "M")],
    estimator=estimator_M_known, est.args = list(q=q),
    generator=generator_M_known, gen.args = list(n=n, p=p, q=q, nmiss=nmiss, missType=missType),
    theta.init <- theta.init,
    H = 20,
    step = 1,
    maxit = 20,
    save=F)
  theta_to_par(fit$theta, p, q)$mu
}

# bias-corrected, unknown M
est_M_unknown <- function(dat){
  par.init <- dat$args[c("A", "psi", "mu")]
  theta.init <- get_theta_init(par.init)
  fit <- ib(dat,
    estimator=estimator_M_unknown,
    est.args=list(q=q),
    generator=generator_M_unknown,
    gen.args=list(n=n, p=p*2, q=q),
    theta.init = theta.init,
    H=20,
    step=.5,
    maxit=20,
    save=F)
  theta_to_par(fit$theta, 2*p, q)$mu[1:3]
}



rep <- 500
nmiss <- .5 * n
missType <- "MCAR"
sim_MCAR <- t(sapply(1:rep, function(s){
  dat <- gen_dat(13123+s)
  c(est_naive_impute(dat),
  est_M_known(dat),
  est_M_unknown(dat))
}))
save(sim_MCAR, file="R/sims/sim_MCAR.Rdata")


nmiss <- .5 * n
missType <- "trunc"
sim_trunc <- t(sapply(1:rep, function(s){
  dat <- gen_dat(13123+s)
  c(est_naive_impute(dat),
  est_M_known(dat),
  est_M_unknown(dat))
}))
save(sim_trunc, file="R/sims/sim_trunc.Rdata")
