library(rstan)
library(arm)
library(parallel)
library(ape)
library(stringr)
library(mvnfast)
library(plyr)
library(reshape2)

RNGkind(kind = "L'Ecuyer-CMRG")
seed <- 420
nsim <- 1000


# for poisson model
deg.coef <- llply(degfit, function(x) extract(x, permuted = TRUE))
mu <- list()
for(ii in 1:nsim) {
  n <- size[[1]]
  inc <- sample(deg.coef$beta_inter, 1)
  sz <- sample(deg.coef$beta_mass, 1)
  mv <- deg.coef$beta_move[sample(size[[1]], 1), ]
  di <- deg.coef$beta_diet[sample(size[[1]], 1), ]

  oo <- c()
  for(jj in seq(size[[1]])) {
    reg <- inc + sz * ecol[[1]]$mass[jj] + 
           sum(mv * move[[1]][jj, ]) + sum(di * diet[[1]][jj, ])
    oo[jj] <- rpois(1, lambda = exp(reg))
  }
  mu[[ii]] <- oo
}

par(mfrow = c(5, 4), mar = c(4, 4, 2, 2))
hist(deg[[1]])
for(s in 1:19)
  hist(mu[[s]])



# for neg binom model
over.coef <- llply(overfit, function(x) extract(x, permuted = TRUE))

#rnbinom(n, size = phi, mu = mu)

