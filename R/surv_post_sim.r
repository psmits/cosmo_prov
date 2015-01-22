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

load('../data/survival_out.rdata')

wei.surv <- function(time, scale, shape) {
  exp(-((1 / scale) * time)**shape)
}

cum.haz <- function(time, scale, shape) {
  -log(wei.surv(time, scale, shape))
}

martingale.res <- function(time, scale, shape, inclusion) {
  inclusion - cum.haz(time, scale, shape)
}

deviance.res <- function(time, scale, shape, inclusion) {
  martin <- martingale.res(time, scale, shape, inclusion)
  si <- sign(martin)
  inn <- martin + (inclusion * log(inclusion - martin))
  under <- -2 * inn
  out <- si * sqrt(under)
  out
}


# zero posterior data set simulations
set.seed(seed)
zpost <- extract(zfit, permuted = TRUE)
zz <- list()
for(i in 1:nsim) {
  oo <- rweibull(length(duration), shape = sample(zpost$alpha, 1), 
                 scale = sample(zpost$sigma, 1))
  zz[[i]] <- oo
}


# weibull posterior data set simulations
mpost <- extract(mfit, permuted = TRUE)
scale.mpost <- extract(scale.mfit, permuted = TRUE)

dat <- cbind(c(scale.data$occ_unc, scale.data$occ_cen), 
             c(scale.data$size_unc, scale.data$size_cen))
dd <- rbind(scale.data$diet_unc, scale.data$diet_cen)
mo <- rbind(scale.data$move_unc, scale.data$move_cen)

set.seed(seed)
dead <- duration
mm <- list()
mm.res <- list()
for(i in 1:nsim) {
  n <- length(dead)
  alp <- sample(scale.mpost$alpha, 1)
  inc <- sample(scale.mpost$beta_inter, 1)
  oc <- sample(scale.mpost$beta_occ, 1)
  sz <- sample(scale.mpost$beta_size, 1)
  mv <- scale.mpost$beta_move[sample(nrow(scale.mpost$beta_move), 1), ]
  di <- scale.mpost$beta_diet[sample(nrow(scale.mpost$beta_diet), 1), ]
  ff <- scale.mpost$rando[sample(nrow(scale.mpost$rando), 1), ]

  oo <- c()
  rr <- c()
  for(j in seq(n)) {
    reg <- inc + oc * dat[j, 1] + sz * dat[j, 2] + 
    sum(di * dd[j, ]) + sum(mv * mo[j, ]) + ff[coh[j]]
    oo[j] <- rweibull(1, scale = exp(-(reg) / alp),
                      shape = alp)
    
    rr[j] <- deviance.res(duration[j], 
                          scale = exp(-(reg) / alp), 
                          shape = alp, 
                          inclusion = extinct[j])
  }
  mm[[i]] <- oo

  mm.res[[i]] <- rr
}


# exponential data set simulations
set.seed(seed)
epost <- extract(efit, permuted = TRUE)
dat <- cbind(c(data$occ_unc, data$occ_cen), 
             c(data$size_unc, data$size_cen))
dd <- rbind(data$diet_unc, data$diet_cen)
mo <- rbind(data$move_unc, data$move_cen)

dead <- duration
ee <- list()
ee.res <- list()
for(i in 1:nsim) {
  n <- length(dead)
  inc <- sample(epost$beta_inter, 1)
  oc <- sample(epost$beta_occ, 1)
  sz <- sample(epost$beta_size, 1)
  mv <- epost$beta_move[sample(nrow(epost$beta_move), 1), ]
  di <- epost$beta_diet[sample(nrow(epost$beta_diet), 1), ]
  ff <- epost$rando[sample(nrow(epost$rando), 1), ]

  oo <- c()
  rr <- c()
  for(j in seq(n)) {
    reg <- inc + oc * dat[j, 1] + sz * dat[j, 2] + 
    sum(di * dd[j, ]) + sum(mv * mo[j, ]) + ff[coh[j]]
    oo[j] <- rexp(1, rate = exp(reg))
    
    rr[j] <- deviance.res(duration[j], 
                          scale = 1 / exp(reg), 
                          shape = 1, 
                          inclusion = extinct[j])
  }
  ee[[i]] <- oo

  ee.res[[i]] <- rr
}


# weibull + phylogenetic scaled data set simulations
set.seed(seed)
phypost <- extract(phy.scalemfit, permuted = TRUE)

dat <- cbind(c(scale.data$occ_unc, scale.data$occ_cen), 
             c(scale.data$size_unc, scale.data$size_cen))
dd <- rbind(scale.data$diet_unc, scale.data$diet_cen)
mo <- rbind(scale.data$move_unc, scale.data$move_cen)

dead <- duration
pm <- list()
pm.res <- list()
for(i in 1:nsim) {
  n <- length(dead)
  alp <- sample(phypost$alpha, 1)
  inc <- sample(phypost$beta_inter, 1)
  oc <- sample(phypost$beta_occ, 1)
  sz <- sample(phypost$beta_size, 1)
  mv <- phypost$beta_move[sample(nrow(phypost$beta_move), 1), ]
  di <- phypost$beta_diet[sample(nrow(phypost$beta_diet), 1), ]
  ff <- phypost$rando[sample(nrow(phypost$rando), 1), ]
  pp <- phypost$phy[sample(nrow(phypost$phy), 1), ]

  oo <- c()
  rr <- c()
  for(j in seq(n)) {
    reg <- inc + oc * dat[j, 1] + sz * dat[j, 2] + 
           sum(di * dd[j, ]) + sum(mv * mo[j, ]) + 
           ff[coh[j]] + pp[j]
    oo[j] <- rweibull(1, scale = exp(-(reg) / alp),
                      shape = alp)
    rr[j] <- deviance.res(duration[j], 
                          scale = exp(-(reg) / alp), 
                          shape = alp, 
                          inclusion = extinct[j])

  }
  pm[[i]] <- oo

  pm.res[[i]] <- rr
}


wei.var <- function(scale, shape) {
  scale**2 * (gamma(1 + (2 / shape)) - (gamma(1 + (1 / shape)))**2)
}
sim.var <- function(n) {
  c.star <- rnorm(n, 0, sample(phypost$fv, 1))
#  p.star <- rnorm(n, 0, sqrt(sample(phypost$sq_sigma, 1) *
#                             sample(diag(tree.vcv), 1)))
  p.star <- rnorm(n, 0, sample(phypost$sigma_phy, 1))

  aa <- sample(phypost$alpha, 1)
  inter <- sample(phypost$beta_inter, 1)
  ss <- exp(-(inter + c.star) / aa)
  yy <- exp(-(inter + p.star) / aa)

  out <- c(y = mean(c(mean(wei.var(ss, aa)), mean(wei.var(yy, aa)))),
           p = var(yy),
           c = var(ss))
  out
}
var.star <- mclapply(1:nsim, function(x) sim.var(50000), 
                     mc.cores = detectCores())

var.star <- data.frame(Reduce(rbind, var.star))
