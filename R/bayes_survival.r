library(rstan)
library(boot)
library(parallel)

source('../R/surv_setup.r')

RNGkind(kind = "L'Ecuyer-CMRG")
seed <- 420

# compile model
zero.weibull <- stan(file = '../stan/zero_weibull.stan')
weibull.model <- stan(file = '../stan/weibull_survival.stan')

# set up variables
duration <- dur
extinct <- ext

# continuous variables
size <- na.ecol$mass
occ <- occ

data <- list(duration = duration,
             size = log(size),
             occ = log(occ),
             diet = factor(na.ecol$diet),
             move = factor(na.ecol$move),
             extinct = extinct)

# break up based on censored or not
grab <- data$extinct == 1
unc <- lapply(data, function(x) x[grab])
cen <- lapply(data, function(x) x[!grab])

data <- list(dur_unc = unc$duration,
             size_unc = unc$size,
             occ_unc = unc$occ,
             diet_unc = unc$diet,
             move_unc = unc$move,
             N_unc = length(unc$duration),
             dur_cen = cen$duration,
             size_cen = cen$size,
             occ_cen = cen$occ,
             diet_cen = cen$diet,
             move_cen = cen$move,
             N_cen = length(cen$duration))

data$diet_unc <- model.matrix( ~ data$diet_unc - 1)[, -1]
data$diet_cen <- model.matrix( ~ data$diet_cen - 1)[, -1]
data$move_unc <- model.matrix( ~ data$move_unc - 1)[, -1]
data$move_cen <- model.matrix( ~ data$move_cen - 1)[, -1]

data$D <- length(unique(na.ecol$diet)) - 1
data$M <- length(unique(na.ecol$move)) - 1
data$N <- data$N_unc + data$N_cen
data$L <- min(c(data$dur_unc, data$dur_cen))
data$samp_unc <- seq(data$N_unc)
data$samp_cen <- seq(from = data$N_unc + 1, 
                     to = data$N_unc + data$N_cen, 
                     by = 1)

exdat <- data[c('N_unc', 'N_cen', 'dur_unc', 'dur_cen', 'L')]

# zero model 
zerolist <- mclapply(1:4, mc.cores = detectCores(),
                     function(x) stan(fit = zero.weibull, 
                                      seed = seed,
                                      data = data,
                                      chains = 1, chain_id = x,
                                      refresh = -1))

zfit <- sflist2stanfit(zerolist)

## larger
#modlist <- mclapply(1:4, mc.cores = detectCores(),
#                     function(x) stan(fit = weibull.model, 
#                                      seed = seed,
#                                      data = data,
#                                      chains = 1, chain_id = x,
#                                      refresh = -1))
#
#mfit <- sflist2stanfit(modlist)
