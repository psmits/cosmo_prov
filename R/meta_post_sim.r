library(rstan)
library(arm)
library(parallel)
library(ape)
library(stringr)
library(plyr)
library(reshape2)
library(igraph)

RNGkind(kind = "L'Ecuyer-CMRG")
seed <- 420
nsim <- 1000

source('../R/metacommunity.r')

# begin by reading each of the pois samples
np <- list.files('../data/meta_dump/meta_sample/', pattern = 'pois_samp')
#apr <- str_extract(np, '\\d+')
#np <- np[order(as.numeric(apr))]
nm <- length(np) / 4
pois.post <- list()
for(ii in seq(nm)) {
  pat <- paste0('pois_samp_', ii, '_[0-9]')
  pois.fit <- list.files('../data/meta_dump/meta_sample', 
                         pattern = pat, full.names = TRUE)
  pois.post[[ii]] <- read_stan_csv(pois.fit)
}

np <- list.files('../data/meta_dump/meta_sample/', pattern = 'nb_samp')
nm <- length(np) / 4
nb.post <- list()
for(ii in seq(nm)) {
  pat <- paste0('nb_samp_', ii, '_[0-9]')
  nb.fit <- list.files('../data/meta_dump/meta_sample', 
                         pattern = pat, full.names = TRUE)
  nb.post[[ii]] <- read_stan_csv(nb.fit)
}


# get some WAIC values
# rough estimate of comparative goodness of fit
colVars <- function (a){
  diff <- a - matrix(colMeans(a), nrow(a), ncol(a), byrow = TRUE)
  vars <- colMeans(diff^2) * nrow(a) / (nrow(a) - 1)
  return(vars)
}
waic <- function (stanfit) {
  log_lik <- extract(stanfit, "log_lik")$log_lik
  lppd <- sum(log(colMeans(exp(log_lik))))
  p_waic_1 <- 2 * sum(log(colMeans(exp(log_lik))) - colMeans(log_lik))
  p_waic_2 <- sum(colVars(log_lik))
  waic_2 <- -2 * lppd + 2 * p_waic_2
  return (list(waic = waic_2, p_waic = p_waic_2, 
               lppd = lppd, p_waic_1 = p_waic_1))
}
pois.waic <- llply(pois.post, waic)
negbin.waic <- llply(nb.post, waic)
waic.diff <- Reduce(rbind, Map(function(x, y) x$waic - y$waic, 
                               pois.waic, negbin.waic))


# posterior simulations
emp <- llply(data, function(x) quantile(x$degree, c(0.25, 0.5, 0.75)))

# poisson model
poisson.sim <- function(fit, datas, nsim) {
  fit <- deg.coef[[9]]
  datas <- data[[9]]
  nsim <- nsim
  
  mu <- list()
  for(ii in 1:nsim) {
    n <- datas$N
    inc <- sample(fit$beta_inter, 1)
    sz <- sample(fit$beta_mass, 1)
    mv <- fit$beta_move[sample(n, 1), ]
    di <- fit$beta_diet[sample(n, 1), ]
    ph <- fit$phy[sample(nrow(fit$phy), 1), ]

    oo <- c()
    for(jj in seq(n)) {
      reg <- inc + sz * datas$mass[jj] + 
      sum(mv * datas$move[jj, ]) + sum(di * datas$diet[jj, ]) + 
      ph[jj] + log(datas$off[jj])
      oo[jj] <- rpois(1, lambda = exp(reg))
    }
    mu[[ii]] <- oo
  }
  mu
}
# this section needs to be heavily modified to deal with the sampling files
deg.coef <- llply(pois.post, function(x) extract(x, permuted = TRUE))
pois.out <- Map(function(x, y) poisson.sim(x, y, nsim = nsim),
                x = deg.coef, y = data)
pois.qout <- llply(pois.out, function(y)
                   laply(y, function(x) 
                         quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)))
pois.cout <- Map(function(x, y) aaply(x, .margins = 1, function(a) a - y),
                 x = pois.qout, y = emp)

ll <- length(deg.coef)
pois.for <- Map(function(x, y) poisson.sim(x, y, nsim = nsim),
                x = deg.coef[-1], y = data[1:(ll - 1)])
pois.qfor <- llply(pois.for, function(y)
                   laply(y, function(x) 
                         quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)))
pois.cfor <- Map(function(x, y) aaply(x, .margins = 1, function(a) a - y),
                 x = pois.qfor, y = emp[1:(test - 1)])


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
    ph <- fit$phy[sample(nrow(fit$phy), 1), ]

    oo <- c()
    for(jj in seq(n)) {
      reg <- inc + sz * data$mass[jj] + 
             sum(mv * data$move[jj, ]) + sum(di * data$diet[jj, ]) +
             ph[jj] + log(data$off[jj])
      oo[jj] <- rnbinom(1, mu = exp(reg), size = phi)
    }
    mu[[ii]] <- oo
  }
  mu
}
over.coef <- llply(nb.post, function(x) extract(x, permuted = TRUE))
negbinom.out <- Map(function(x, y) negbinom.sim(x, y, nsim = nsim),
                    x = over.coef, y = data)
negbinom.qout <- llply(negbinom.out, function(y)
                       laply(y, function(x) 
                             quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)))
negbinom.cout <- Map(function(x, y) aaply(x, .margins = 1, function(a) a - y),
                     x = negbinom.qout, y = emp)

ll <- length(over.coef)
negbinom.for <- Map(function(x, y) negbinom.sim(x, y, nsim = nsim),
                    x = over.coef[-1], y = data[1:(ll - 1)])
negbinom.qfor <- llply(negbinom.for, function(y)
                       laply(y, function(x) 
                             quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)))
negbinom.cfor <- Map(function(x, y) aaply(x, .margins = 1, function(a) a - y),
                     x = negbinom.qfor, y = emp[1:(ll - 1)])
