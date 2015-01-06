library(rstan)
library(arm)
library(parallel)
library(ape)
library(stringr)
library(plyr)
library(reshape2)

RNGkind(kind = "L'Ecuyer-CMRG")
seed <- 420
nsim <- 1000


# get some WAIC values
# rough estimate of comparative goodness of fit
colVars <- function (a){
  diff <- a - matrix (colMeans(a), nrow(a), ncol(a), byrow=TRUE)
  vars <- colMeans (diff^2)*nrow(a)/(nrow(a)-1)
  return (vars)
}
waic <- function (stanfit){
  log_lik <- extract (stanfit, "log_lik")$log_lik
  lppd <- sum (log (colMeans(exp(log_lik))))
  p_waic_1 <- 2*sum (log(colMeans(exp(log_lik))) - colMeans(log_lik))
  p_waic_2 <- sum (colVars(log_lik))
  waic_2 <- -2*lppd + 2*p_waic_2
  return (list (waic=waic_2, p_waic=p_waic_2, lppd=lppd, p_waic_1=p_waic_1))
}
pois.waic <- llply(degfit, waic)
negbin.waic <- llply(overfit, waic)


# posterior simulations
# poisson model
poisson.sim <- function(fit, data, nsim) {
  mu <- list()
  for(ii in 1:nsim) {
    n <- data$N
    inc <- sample(fit$beta_inter, 1)
    sz <- sample(fit$beta_mass, 1)
    mv <- fit$beta_move[sample(n, 1), ]
    di <- fit$beta_diet[sample(n, 1), ]

    oo <- c()
    for(jj in seq(n)) {
      reg <- inc + sz * data$mass[jj] + 
      sum(mv * data$move[jj, ]) + sum(di * data$diet[jj, ])
      oo[jj] <- rpois(1, lambda = exp(reg))
    }
    mu[[ii]] <- oo
  }
  mu
}

# remember, the order is reversed
# 2, 4, 6, etc.
deg.coef <- llply(degfit, function(x) extract(x, permuted = TRUE))
pois.out <- Map(function(x, y) poisson.sim(x, y, nsim = nsim),
                x = deg.coef, y = data)

ll <- length(deg.coef)
pois.for <- Map(function(x, y) poisson.sim(x, y, nsim = nsim),
                x = deg.coef[-1], y = data[-ll])

ps.forward <- Map(function(x, y) quantile(x, 0.5) - quantile(y, 0.5),
                  x = pois.out, y = pois.for)


# neg binom model
negbinom.sim <- function(fit, data, nsim) {
  mu <- list()
  for(ii in 1:nsim) {
    n <- data$N
    inc <- sample(fit$beta_inter, 1)
    sz <- sample(fit$beta_mass, 1)
    mv <- fit$beta_move[sample(n, 1), ]
    di <- fit$beta_diet[sample(n, 1), ]
    phi <- sample(fit$phi, 1)

    oo <- c()
    for(jj in seq(n)) {
      reg <- inc + sz * data$mass[jj] + 
      sum(mv * data$move[jj, ]) + sum(di * data$diet[jj, ])
      oo[jj] <- rnbinom(1, mu = exp(reg), size = phi)
    }
    mu[[ii]] <- oo
  }
  mu
}
over.coef <- llply(overfit, function(x) extract(x, permuted = TRUE))
negbinom.out <- Map(function(x, y) negbinom.sim(x, y, nsim = nsim),
                    x = over.coef, y = data)

ll <- length(over.coef)
negbinom.for <- Map(function(x, y) negbinom.sim(x, y, nsim = nsim),
                    x = over.coef[-1], y = data[-ll])

nb.forward <- Map(function(x, y) quantile(x, 0.5) - quantile(y, 0.5),
                  x = negbinom.out, y = negbinom.for)
