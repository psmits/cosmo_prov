library(ggplot2)
library(reshape2)
library(scales)
library(hexbin)
library(stringr)
library(grid)

theme_set(theme_bw())
cbp <- c('#E69F00', '#56B4E9', '#009E73', '#F0E442', 
         '#0072B2', '#D55E00', '#CC79A7')
theme_update(axis.text = element_text(size = 20),
             axis.title = element_text(size = 30),
             axis.title.y = element_text(hjust = 0.1),
             legend.text = element_text(size = 25),
             legend.title = element_text(size = 26),
             legend.key.size = unit(2, 'cm'),
             strip.text = element_text(size = 25))

# marginal posteriors
melted <- melt(mpost)
ns <- length(mpost$lp__)

# create the different combinations
regcoef <- melted[str_detect(melted$L1, 'beta'), ]
base.inter <- melted[melted$L1 == 'beta_inter', ]  # base line
move.eff <- melted[melted$L1 == 'beta_move', ]  # move effect
diet.eff <- melted[melted$L1 == 'beta_diet', ]  # diet effect

# make fancy versions of these
# lower values, lower risk
hist(c(base.inter$value + move.eff$value)[1:ns] - base.inter$value) # ground - arb
hist(c(base.inter$value + move.eff$value)[(ns+1):(2*ns)] - base.inter$value) # scan - arb 
hist(c(base.inter$value + move.eff$value)[(ns+1):(2*ns)] - 
     c(base.inter$value + move.eff$value[1:ns])) # scan - ground

hist(c(base.inter$value + diet.eff$value)[1:ns] - 
     base.inter$value)  # herb - carn
hist(c(base.inter$value + diet.eff$value)[(ns+1):(2*ns)] - 
     base.inter$value) # insect - carn
hist(c(base.inter$value + diet.eff$value)[(2*ns+1):(3*ns)] - 
     base.inter$value) # omni - carn
hist(c(base.inter$value + diet.eff$value)[(2*ns+1):(3*ns)] - 
     c(base.inter$value + diet.eff$value)[1:ns]) # omni - herb
hist(c(base.inter$value + diet.eff$value)[(2*ns+1):(3*ns)] - 
     c(base.inter$value + diet.eff$value)[(ns+1):(2*ns)]) # omni - insect



# histogram of variance in cohort effect
# for cohort effect, do point range with 80% quartile
