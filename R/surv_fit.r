library(rstan)
library(ggplot2)
library(scales)
library(grid)
library(arm)
library(parallel)
library(xtable)
library(ape)
library(stringr)

janis <- read.csv('../data/suspicious.csv', header = FALSE)

RNGkind(kind = "L'Ecuyer-CMRG")
seed <- 420
nsim <- 1000

set.seed(seed)
source('../R/surv_setup.r')  # data cleaned up nice

janis.ecol <- na.ecol[na.ecol$taxa %in% janis[, 1], ]
janis.dur <- dur[na.ecol$taxa %in% janis[, 1]]
janis.ext <- ext[na.ecol$taxa %in% janis[, 1]]
janis.coh <- cohort[na.ecol$taxa %in% janis[, 1]]
janis.occ <- occ[na.ecol$taxa %in% janis[, 1]]

nojanis.ecol <- na.ecol[!(na.ecol$taxa %in% janis[, 1]), ]
nojanis.dur <- dur[!(na.ecol$taxa %in% janis[, 1])]
nojanis.ext <- ext[!(na.ecol$taxa %in% janis[, 1])]
nojanis.coh <- cohort[!(na.ecol$taxa %in% janis[, 1])]
nojanis.occ <- occ[!(na.ecol$taxa %in% janis[, 1])]


# super tree
load('../data/scaled_super.rdata')  # best

plio <- which(cohort == 1)

# set up variables
duration <- dur[-plio]
extinct <- ext[-plio]
coh <- cohort[-plio] - 1

nojanis.duration <- nojanis.dur[-plio]
nojanis.extinct <- nojanis.ext[-plio]
nojanis.coh <- nojanis.coh[-plio] 

# continuous variables
size <- na.ecol$mass[-plio]
occ <- occ[-plio]
di <- factor(na.ecol$diet[-plio])
mo <- factor(na.ecol$move[-plio])
inr <- interaction(di, mo)

nojanis.size <- nojanis.ecol$mass[-plio]
nojanis.occ <- nojanis.occ[-plio]
nojanis.di <- factor(nojanis.ecol$diet[-plio])
nojanis.mo <- factor(nojanis.ecol$move[-plio])
nojanis.inr <- interaction(nojanis.di, nojanis.mo)


# compare distribution of body sizes
theme_set(theme_bw())
cbp <- c('#E69F00', '#56B4E9', '#009E73', '#F0E442', 
         '#0072B2', '#D55E00', '#CC79A7')
theme_update(axis.text = element_text(size = 6),
             axis.title = element_text(size = 10),
             legend.text = element_text(size = 8),
             legend.title = element_text(size = 9),
             legend.key.size = unit(0.75, 'cm'),
             strip.text = element_text(size = 8))
jj <- na.ecol[-plio, 1] %in% as.character(janis[, 1])
bs <- data.frame(size = log(size),
                 state = ifelse(jj == 1, 'suspicious', 'fine'),
                 time = coh * 2,
                 diet = di,
                 move = mo)
bs.gg <- ggplot(bs, aes(x = size, fill = state))
bs.gg <- bs.gg + geom_histogram()
bs.gg <- bs.gg + scale_fill_manual(values = c('black', 'blue'))
bs.gg <- bs.gg + labs(x = 'log(body size in grams)', y = 'Frequency')
ggsave(filename = '../doc/na_surv/figure/body_size_compare.png', plot = bs.gg,
       width = 4, height = 3)
bs.gg <- bs.gg + facet_wrap( ~ time)
ggsave(filename = '../doc/na_surv/figure/body_size_compare_time.png', plot = bs.gg,
       width = 8, height = 6)

# prepare the vcv
tax.nam <- str_replace(na.ecol[-plio, 1], ' ', '_')
to.drop <- spt$tip.label[!(spt$tip.label %in% tax.nam)]
na.tree <- drop.tip(spt, to.drop)
na.tree$edge.length <- na.tree$edge.length / max(diag(vcv(na.tree)))
tree.vcv <- vcv(na.tree)
split.tax <- rev(split(tax.nam, extinct))
cor.ord <- match(unlist(split.tax), colnames(tree.vcv))
tree.vcv <- tree.vcv[cor.ord, cor.ord]


# prepare the _other_ vcv
nj.tax.nam <- str_replace(nojanis.ecol[-plio, 1], ' ', '_')
nj.to.drop <- spt$tip.label[!(spt$tip.label %in% nj.tax.nam)]
nj.na.tree <- drop.tip(spt, nj.to.drop)
nj.na.tree$edge.length <- nj.na.tree$edge.length / max(diag(vcv(nj.na.tree)))
nj.tree.vcv <- vcv(nj.na.tree)
nj.split.tax <- rev(split(nj.tax.nam, nojanis.extinct))
nj.cor.ord <- match(unlist(nj.split.tax), colnames(nj.tree.vcv))
nj.tree.vcv <- tree.vcv[nj.cor.ord, nj.cor.ord]


#length(colnames(tree.vcv))
#sum(str_replace(colnames(tree.vcv), '_', ' ') %in% as.character(north.source[north.source$Source == 'PBDB + regression', 1]))
#cov2cor(tree.vcv)

# full dataset
data <- list(duration = duration,
             size = rescale(log(size)),
             occ = rescale(logit(occ)),
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
with(data, {stan_rdump(list = c('N', 'N_unc', 'N_cen', 'D', 'M', 'L', 'C', 
                                'dur_unc', 'dur_cen', 'size_unc', 'size_cen', 
                                'occ_unc', 'occ_cen', 'diet_unc', 'diet_cen', 
                                'move_unc', 'move_cen', 'coh_unc', 'coh_cen', 
                                'samp_unc', 'samp_cen', 'vcv'),
                       file = '../data/data_dump/surv_info.data.R')})

# w/o suspicious
data <- list(duration = nojanis.duration,
             size = rescale(log(nojanis.size)),
             occ = rescale(logit(nojanis.occ)),
             diet = nojanis.di,
             move = nojanis.mo,
             rac = nojanis.inr,
             extinct = nojanis.extinct,
             coh = nojanis.coh)

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

data$D <- length(unique(nojanis.ecol$diet)) - 1
data$M <- length(unique(nojanis.ecol$move)) - 1
data$R <- length(unique(nojanis.inr)) - 1
data$N <- data$N_unc + data$N_cen
data$L <- min(c(data$dur_unc, data$dur_cen))
data$C <- length(unique(nojanis.coh))
data$samp_unc <- seq(data$N_unc)
data$samp_cen <- seq(from = data$N_unc + 1, 
                     to = data$N_unc + data$N_cen, 
                     by = 1)
data$vcv <- nj.tree.vcv

# rescale following gelman
with(data, {stan_rdump(list = c('N', 'N_unc', 'N_cen', 'D', 'M', 'L', 'C', 
                                'dur_unc', 'dur_cen', 'size_unc', 'size_cen', 
                                'occ_unc', 'occ_cen', 'diet_unc', 'diet_cen', 
                                'move_unc', 'move_cen', 'coh_unc', 'coh_cen', 
                                'samp_unc', 'samp_cen', 'vcv'),
                       file = '../data/data_dump/surv_nojanis_info.data.R')})
