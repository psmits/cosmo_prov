library(rstan)
library(arm)
library(plyr)
library(reshape2)
library(stringr)
library(xtable)

source('../R/surv_post_sim.r')

exct <- extract(phy.scalemfit)
pp <- !(names(exct) %in% c('rando', 'phy', 'sigma_phy', 'lp__', 'log_lik'))
summs <- summary(phy.scalemfit, pars = names(exct)[pp])[[1]]
summs <- summs[, !colnames(summs) %in% c('se_mean', 'n_eff')]

rownames(summs) <- c('intercept', 'logit(occupancy)', 'log(size)', 
                     'alpha', 'ground dwelling', 'scansorial', 
                     'herbivore', 'insectivore', 'omnivore',
                     'sd cohort', 'sd phylogeny')
summs <- summs[c(4, 1:3, 5:nrow(summs)), ]

sum.table <- xtable(summs, label = 'post_sum')
print.xtable(sum.table, file = '../doc/na_surv/post_out_raw_noj.tex')
