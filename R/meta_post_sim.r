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


# takes posteriors and data
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

# for poisson model
deg.coef <- llply(degfit, function(x) extract(x, permuted = TRUE))
pois.out <- Map(function(x, y) poisson.sim(x, y, nsim = nsim),
                x = deg.coef, y = data)
ll <- length(deg.coef)
pois.for <- Map(function(x, y) poisson.sim(x, y, nsim = nsim),
                x = deg.coef[-ll], y = data[-1])

#par(mfrow = c(5, 4), mar = c(4, 4, 2, 2))
#hist(data[[1]]$deg)
#for(s in 1:19)
#  hist(pois.out[[s]])


# neg binom simulations
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
                    x = over.coef[-ll], y = data[-1])
