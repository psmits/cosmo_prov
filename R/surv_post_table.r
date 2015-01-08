library(rstan)
library(arm)
library(plyr)
library(reshape2)
library(stringr)
library(xtable)

exct <- extract(phy.scalemfit)
pp <- !(names(exct) %in% c('rando', 'phy', 'sigma_phy', 'lp__'))
summs <- summary(phy.scalemfit, pars = names(exct)[pp])[[1]]
summs <- summs[rownames(summs) != 'alpha', 
               !colnames(summs) %in% c('se_mean', 'n_eff')]

rownames(summs) <- c('intercept', 'logit(occupancy)', 'log(size)', 
                     'ground dwelling', 'scansorial', 
                     'herbivore', 'insectivore', 'omnivore',
                     'sd cohort', 'sd phylogeny')

sum.table <- xtable(summs, label = 'post_sum')
print.xtable(sum.table, file = '../doc/na_surv/post_out_raw.tex')
