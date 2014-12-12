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
zfit.table <- xtable(summary(zfit)[[1]], label = 'zfit_tab')
print.xtable(zfit.table, 
             file = '../doc/na_surv/zfit_table_raw.tex')

set.seed(seed)
zpost <- extract(zfit, permuted = TRUE)
zz <- list()
for(i in 1:nsim) {
  oo <- rweibull(length(duration), shape = sample(zpost$alpha, 1), 
                 scale = sample(zpost$sigma, 1))
  zz[[i]] <- oo
}


# larger
modlist <- mclapply(1:4, mc.cores = detectCores(),
                    function(x) stan(fit = weibull.model, 
                                     seed = seed,
                                     data = data,
                                     chains = 1, chain_id = x,
                                     refresh = -1))

mfit <- sflist2stanfit(modlist)
mfit.table <- xtable(summary(mfit)[[1]], label = 'weifit_tab')
print.xtable(mfit.table, 
             file = '../doc/na_surv/mfit_table_raw.tex')


scale.modlist <- mclapply(1:4, mc.cores = detectCores(),
                          function(x) stan(fit = weibull.model, 
                                           seed = seed,
                                           data = scale.data,
                                           chains = 1, chain_id = x,
                                           refresh = -1))

scale.mfit <- sflist2stanfit(scale.modlist)
scalefit.table <- xtable(summary(scale.mfit)[[1]], label = 'scalefit_tab')
print.xtable(scalefit.table, 
             file = '../doc/na_surv/scalefit_table_raw.tex')


explist <- mclapply(1:4, mc.cores = detectCores(),
                    function(x) stan(fit = exponential.model, 
                                     seed = seed,
                                     data = data,
                                     chains = 1, chain_id = x,
                                     refresh = -1))

efit <- sflist2stanfit(explist)
efit.table <- xtable(summary(efit)[[1]], label = 'efit_tab')
print.xtable(efit.table, 
             file = '../doc/na_surv/efit_table_raw.tex')


mpost <- extract(mfit, permuted = TRUE)
scale.mpost <- extract(scale.mfit, permuted = TRUE)

dat <- cbind(c(data$occ_unc, data$occ_cen), 
             c(data$size_unc, data$size_cen))
dd <- rbind(data$diet_unc, data$diet_cen)
mo <- rbind(data$move_unc, data$move_cen)

set.seed(seed)
dead <- duration
mm <- list()
for(i in 1:nsim) {
  n <- length(dead)
  alp <- sample(mpost$alpha, 1)
  inc <- sample(mpost$beta_inter, 1)
  oc <- sample(mpost$beta_occ, 1)
  sz <- sample(mpost$beta_size, 1)
  mv <- mpost$beta_move[sample(nrow(mpost$beta_move), 1), ]
  di <- mpost$beta_diet[sample(nrow(mpost$beta_diet), 1), ]
  ff <- mpost$rando[sample(nrow(mpost$rando), 1), ]

  oo <- c()
  for(j in seq(n)) {
    reg <- inc + oc * dat[j, 1] + sz * dat[j, 2] + 
    sum(di * dd[j, ]) + sum(mv * mo[j, ]) + ff[coh[j]]
    oo[j] <- rweibull(1, scale = exp(-reg) / alp,
                      shape = alp)
  }
  mm[[i]] <- oo
}


set.seed(seed)
epost <- extract(efit, permuted = TRUE)
dat <- cbind(c(data$occ_unc, data$occ_cen), 
             c(data$size_unc, data$size_cen))
dd <- rbind(data$diet_unc, data$diet_cen)
mo <- rbind(data$move_unc, data$move_cen)

dead <- duration
ee <- list()
for(i in 1:nsim) {
  n <- length(dead)
  inc <- sample(epost$beta_inter, 1)
  oc <- sample(epost$beta_occ, 1)
  sz <- sample(epost$beta_size, 1)
  mv <- epost$beta_move[sample(nrow(epost$beta_move), 1), ]
  di <- epost$beta_diet[sample(nrow(epost$beta_diet), 1), ]
  ff <- epost$rando[sample(nrow(epost$rando), 1), ]

  oo <- c()
  for(j in seq(n)) {
    reg <- inc + oc * dat[j, 1] + sz * dat[j, 2] + 
    sum(di * dd[j, ]) + sum(mv * mo[j, ]) + ff[coh[j]]
    oo[j] <- rexp(1, rate = exp(reg))
  }
  ee[[i]] <- oo
}


# phylogenetic random effect models
#phylist <- mclapply(1:4, mc.cores = detectCores(),
#                    function(x) stan(fit = phywei.model, 
#                                     seed = seed,
#                                     data = data,
#                                     iter = 4000,
#                                     chains = 1, chain_id = x,
#                                     refresh = -1))
#phy.mfit <- sflist2stanfit(phylist)
#
#scale.phylist <- mclapply(1:4, mc.cores = detectCores(),
#                    function(x) stan(fit = phywei.model, 
#                                     seed = seed,
#                                     data = scale.data,
#                                     iter = 4000,
#                                     chains = 1, chain_id = x,
#                                     refresh = -1))
#phy.scalemfit <- sflist2stanfit(scale.phylist)
#
#set.seed(seed)
#phypost <- extract(phy.scalemfit, permuted = TRUE)
#
#dat <- cbind(c(data$occ_unc, data$occ_cen), 
#             c(data$size_unc, data$size_cen))
#dd <- rbind(data$diet_unc, data$diet_cen)
#mo <- rbind(data$move_unc, data$move_cen)
#
#dead <- duration
#pm <- list()
#for(i in 1:100) {
#  n <- length(dead)
#  alp <- sample(phypost$alpha, 1)
#  inc <- sample(phypost$beta_inter, 1)
#  oc <- sample(phypost$beta_occ, 1)
#  sz <- sample(phypost$beta_size, 1)
#  mv <- phypost$beta_move[sample(nrow(phypost$beta_move), 1), ]
#  di <- phypost$beta_diet[sample(nrow(phypost$beta_diet), 1), ]
#  ff <- phypost$rando[sample(nrow(phypost$rando), 1), ]
#  pp <- phypost$phy
#
#  oo <- c()
#  for(j in seq(n)) {
#    reg <- inc + oc * dat[j, 1] + sz * dat[j, 2] + 
#    sum(di * dd[j, ]) + sum(mv * mo[j, ]) + ff[coh[j]] + pp[j]
#    oo[j] <- rweibull(1, scale = exp(-reg) / alp,
#                      shape = alp)
#  }
#  pm[[i]] <- oo
#}
