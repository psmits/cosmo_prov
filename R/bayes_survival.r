library(rstan)
library(boot)
library(parallel)
library(truncdist)

source('../R/surv_setup.r')

RNGkind(kind = "L'Ecuyer-CMRG")
seed <- 420

# compile model
zero.weibull <- stan(file = '../stan/zero_weibull.stan')
weibull.model <- stan(file = '../stan/weibull_survival.stan')

plio <- which(cohort == 1)

# set up variables
duration <- dur[-plio]
extinct <- ext[-plio]
coh <- cohort[-plio] - 1

# continuous variables
size <- na.ecol$mass[-plio]
occ <- occ[-plio]
di <- factor(na.ecol$diet[-plio])
mo <- factor(na.ecol$move[-plio])
inr <- interaction(di, mo)

data <- list(duration = duration,
             size = log(size),
             occ = log(occ),
             diet = di,
             move = mo,
             rac = inr,
             extinct = extinct,
             coh = coh)

# break up based on censored or not
grab <- data$extinct == 1
unc <- lapply(data, function(x) x[grab])
cen <- lapply(data, function(x) x[!grab])

data <- list(dur_unc = unc$duration,
             size_unc = unc$size,
             occ_unc = unc$occ,
             diet_unc = unc$diet,
             move_unc = unc$move,
             coh_unc = unc$coh,
             rac_unc = unc$rac,
             N_unc = length(unc$duration),
             dur_cen = cen$duration,
             size_cen = cen$size,
             occ_cen = cen$occ,
             diet_cen = cen$diet,
             move_cen = cen$move,
             coh_cen = cen$coh,
             rac_cen = cen$rac,
             N_cen = length(cen$duration))

data$diet_unc <- model.matrix( ~ data$diet_unc - 1)[, -1]
data$diet_cen <- model.matrix( ~ data$diet_cen - 1)[, -1]
data$move_unc <- model.matrix( ~ data$move_unc - 1)[, -1]
data$move_cen <- model.matrix( ~ data$move_cen - 1)[, -1]
data$rac_unc <- model.matrix( ~ data$rac_unc - 1)[, -1]
data$rac_cen <- model.matrix( ~ data$rac_cen - 1)[, -1]

data$D <- length(unique(na.ecol$diet)) - 1
data$M <- length(unique(na.ecol$move)) - 1
data$R <- length(unique(inr)) - 1
data$N <- data$N_unc + data$N_cen
data$L <- min(c(data$dur_unc, data$dur_cen))
data$C <- length(unique(coh))
data$samp_unc <- seq(data$N_unc)
data$samp_cen <- seq(from = data$N_unc + 1, 
                     to = data$N_unc + data$N_cen, 
                     by = 1)

exdat <- data[c('N_unc', 'N_cen', 'dur_unc', 'dur_cen', 'L')]

## zero model 
#zerolist <- mclapply(1:4, mc.cores = detectCores(),
#                     function(x) stan(fit = zero.weibull, 
#                                      seed = seed,
#                                      data = data,
#                                      chains = 1, chain_id = x,
#                                      refresh = -1))
#
#zfit <- sflist2stanfit(zerolist)
#
#zpost <- extract(zfit, permuted = TRUE)
#mm <- list()
#for(i in 1:1000) {
#  oo <- rweibull(length(duration), shape = sample(zpost$alpha, 1), 
#                 scale = sample(zpost$sigma, 1))
#  mm[[i]] <- oo
#}
#tp <- sum(laply(mm, mean) > mean(duration))
#br <- seq(0, max(ceiling(laply(mm, max))), 1)
#par(mfrow = c(5, 4), mar = c(4, 4, 2, 2))
#hist(duration, breaks = br)
#for(s in 1:19)
#  hist(mm[[s]], breaks = br)
#
#
# larger
#modlist <- mclapply(1:4, mc.cores = detectCores(),
#                     function(x) stan(fit = weibull.model, 
#                                      seed = seed,
#                                      data = data,
#                                      chains = 1, chain_id = x,
#                                      refresh = -1))
#
#mfit <- sflist2stanfit(modlist)
#
#mpost <- extract(mfit, permuted = TRUE)
#  
#dat <- cbind(c(data$occ_unc, data$occ_cen), 
#             c(data$size_unc, data$size_cen))
#dd <- rbind(data$diet_unc, data$diet_cen)
#mo <- rbind(data$move_unc, data$move_cen)
#
#dead <- duration
#mm <- list()
#for(i in 1:100) {
#  n <- length(dead)
#  alp <- sample(mpost$alpha, 1)
#  inc <- sample(mpost$beta_inter, 1)
#  oc <- sample(mpost$beta_occ, 1)
#  sz <- sample(mpost$beta_size, 1)
#  mv <- mpost$beta_move[sample(nrow(mpost$beta_move), 1), ]
#  di <- mpost$beta_diet[sample(nrow(mpost$beta_diet), 1), ]
#  ff <- mpost$frailty[sample(nrow(mpost$frailty), 1), ]
#
#  oo <- c()
#  for(j in seq(n)) {
#    reg <- inc + oc * dat[j, 1] + sz * dat[j, 2] + 
#           sum(di * dd[j, ]) + sum(mv * mo[j, ]) + ff[coh[j]]
#    oo[j] <- rweibull(1, scale = exp(-reg) / alp,
#                      shape = alp)
#  }
#  mm[[i]] <- oo
#}
#
#sum(laply(mm, function(x) mean(x)) > mean(dead))
#hist(laply(mm, function(x) mean(x))); abline(v = mean(dead))
#sum(laply(mm, function(x) quantile(x, .5)) > quantile(dead, .5))
#hist(laply(mm, function(x) quantile(x, .5))); abline(v = quantile(dead, .5))
#sum(laply(mm, function(x) quantile(x, .75)) > quantile(dead, .75))
#hist(laply(mm, function(x) quantile(x, .75))); abline(v = quantile(dead, .75))
#
#br <- seq(0, max(ceiling(laply(mm, max))), 1)
#par(mfrow = c(5, 4), mar = c(4, 4, 2, 2))
#hist(dead, breaks = br)
#for(s in 1:19)
#  hist(mm[[s]], breaks = br)
