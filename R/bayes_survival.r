library(rstan)
library(boot)
library(parallel)

#source('../R/surv_setup.r')

RNGkind(kind = "L'Ecuyer-CMRG")
seed <- 420

# compile model
weibull.model <- stan(file = '../stan/weibull_survival.stan')

# set up variables
duration <- na.surv[, 1]
extinct <- na.surv[, 2]

# index categories
diets <- model.matrix( ~ diet - 1, data = na.ecol)
herb <- diets[, 2]
insect <- diets[, 3]
omni <- diets[, 4]

moves <- model.matrix( ~ move - 1, data = na.ecol)
ground <- moves[, 2]
scans <- moves[, 3]

# continuous variables
size <- na.ecol$mass
occ <- na.ecol$occ


# break up based on censored or not
data <- list(duration = duration,
             size = log(size),
             occ = occ,
             herb = herb,
             insect = insect,
             omni = omni,
             ground = ground,
             scans = scans,
             extinct = extinct)

grab <- data$extinct == 1
unc <- lapply(data, function(x) x[grab])
cen <- lapply(data, function(x) x[!grab])

data <- list(dur_unc = unc$duration,
             size_unc = unc$size,
             occ_unc = unc$occ,
             herb_unc = unc$herb,
             insect_unc = unc$insect,
             omni_unc = unc$omni,
             ground_unc = unc$ground,
             scans_unc = unc$scans,
             N_unc = length(unc$duration),
             dur_cen = cen$duration,
             size_cen = cen$size,
             occ_cen = cen$occ,
             herb_cen = cen$herb,
             insect_cen = cen$insect,
             omni_cen = cen$omni,
             ground_cen = cen$ground,
             scans_cen = cen$scans,
             N_cen = length(cen$duration))

data$N <- data$N_unc + data$N_cen
data$samp_unc <- seq(data$N_unc)
data$samp_cen <- seq(from = data$N_unc + 1, 
                     to = data$N_unc + data$N_cen, 
                     by = 1)

## weibull
weilist <- mclapply(1:4, mc.cores = detectCores(),
                    function(x) stan(fit = weibull.model, 
                                     seed = seed,
                                     data = data,
                                     chains = 1, chain_id = x,
                                     refresh = -1))

wfit <- sflist2stanfit(weilist)
