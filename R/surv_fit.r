library(rstan)
library(arm)
library(parallel)
library(xtable)
library(ape)
library(stringr)

RNGkind(kind = "L'Ecuyer-CMRG")
seed <- 420
nsim <- 1000

set.seed(seed)
source('../R/surv_setup.r')
load('../data/taxonomy_tree.rdata')

# compile model
zero.weibull <- stan(file = '../stan/zero_weibull.stan')
weibull.model <- stan(file = '../stan/weibull_survival.stan')
phywei.model <- stan(file = '../stan/weibull_phy_surv.stan')
exponential.model <- stan(file = '../stan/exp_survival.stan')

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

# prepare the vcv
tax.nam <- str_replace(na.ecol[-plio, 1], ' ', '_')
to.drop <- na.scale$tip.label[!(na.scale$tip.label %in% tax.nam)]
na.tree <- drop.tip(na.scale, to.drop)
na.tree$edge.length <- na.tree$edge.length / max(diag(vcv(na.tree)))
tree.vcv <- vcv(na.tree)
split.tax <- rev(split(tax.nam, extinct))
cor.ord <- match(unlist(split.tax), colnames(tree.vcv))
tree.vcv <- tree.vcv[cor.ord, cor.ord]
#cov2cor(tree.vcv)

data <- list(duration = duration,
             size = log(size),
             occ = logit(occ),
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
data$vcv <- tree.vcv

# rescale following gelman
scale.data <- data
scale.data$size_unc <- rescale(scale.data$size_unc)
scale.data$size_cen <- rescale(scale.data$size_cen)
scale.data$occ_unc <- rescale(scale.data$occ_unc)
scale.data$occ_cen <- rescale(scale.data$occ_cen)

exdat <- data[c('N_unc', 'N_cen', 'dur_unc', 'dur_cen', 'L')]

# zero model 
zerolist <- mclapply(1:4, mc.cores = detectCores(),
                     function(x) stan(fit = zero.weibull, 
                                      seed = seed,
                                      data = data,
                                      chains = 1, chain_id = x,
                                      refresh = -1))

zfit <- sflist2stanfit(zerolist)


# weibull model minus phylogeny
modlist <- mclapply(1:4, mc.cores = detectCores(),
                    function(x) stan(fit = weibull.model, 
                                     seed = seed,
                                     data = data,
                                     iter = 1000,
                                     chains = 1, chain_id = x,
                                     refresh = -1))

mfit <- sflist2stanfit(modlist)


# same with scaled data
scale.modlist <- mclapply(1:4, mc.cores = detectCores(),
                          function(x) stan(fit = weibull.model, 
                                           seed = seed,
                                           data = scale.data,
                                           iter = 1000,
                                           chains = 1, chain_id = x,
                                           refresh = -1))

scale.mfit <- sflist2stanfit(scale.modlist)


# exponential model minus phylogeny
# for comparison with weibull to see if value of alpha merited
explist <- mclapply(1:4, mc.cores = detectCores(),
                    function(x) stan(fit = exponential.model, 
                                     seed = seed,
                                     data = data,
                                     iter = 1000,
                                     chains = 1, chain_id = x,
                                     refresh = -1))

efit <- sflist2stanfit(explist)


# phylogenetic random effect models
phylist <- mclapply(1:4, mc.cores = detectCores(),
                    function(x) stan(fit = phywei.model, 
                                     seed = seed,
                                     data = data,
                                     iter = 20000,
                                     thin = 20,
                                     chains = 1, chain_id = x,
                                     refresh = -1))
phy.mfit <- sflist2stanfit(phylist)

scale.phylist <- mclapply(1:4, mc.cores = detectCores(),
                          function(x) stan(fit = phywei.model, 
                                           seed = seed,
                                           data = scale.data,
                                           iter = 20000,
                                           thin = 20,
                                           chains = 1, chain_id = x,
                                           refresh = -1))
phy.scalemfit <- sflist2stanfit(scale.phylist)

save.image(file = '../data/survival_out.rdata')

