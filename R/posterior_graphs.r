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
# these are particularilty useful for summary statistics
hist(c(base.inter$value + move.eff$value)[1:ns] - base.inter$value) # ground - arb
hist(c(base.inter$value + move.eff$value)[(ns+1):(2*ns)] - base.inter$value) # scan - arb 
hist(c(base.inter$value + move.eff$value)[(ns+1):(2*ns)] - 
     c(base.inter$value + move.eff$value[1:ns])) # scan - ground

# effect of locomotor category
arb.eff <- base.inter$value
grd.eff <- base.inter$value + move.eff$value[1:ns]
scn.eff <- base.inter$value + move.eff$value[(ns + 1):(2 * ns)]
loceff <- melt(cbind(arb.eff, grd.eff, scn.eff))
loco <- ggplot(loceff, aes(x = value))
loco <- loco + geom_histogram(aes(y = ..density..))
loco <- loco + geom_vline(xintercept = 0, colour = 'grey', size = 2)
loco <- loco + facet_grid(Var2 ~ .)
loco <- loco + labs(x = 'Value', y = 'Prob. Density')


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

# effect of dietary category
crn.eff <- base.inter$value
hrb.eff <- base.inter$value + diet.eff$value[1:ns]
ist.eff <- base.inter$value + diet.eff$value[(ns + 1):(2 * ns)]
omn.eff <- base.inter$value + diet.eff$value[(2 * ns + 1):(3 * ns)]
deteff <- melt(cbind(crn.eff, hrb.eff, ist.eff, omn.eff))
diet <- ggplot(deteff, aes(x = value))
diet <- diet + geom_histogram(aes(y = ..density..))
diet <- diet + geom_vline(xintercept = 0, colour = 'grey', size = 2)
diet <- diet + facet_grid(Var2 ~ .)
diet <- diet + labs(x = 'Value', y = 'Prob. Density')



# histogram of variance in cohort effect
# for cohort effect, do point range with 80% quartile
top <- apply(mpost$rando, 2, function(x) quantile(x, .8))
bot <- apply(mpost$rando, 2, function(x) quantile(x, .2))
med <- apply(mpost$rando, 2, function(x) quantile(x, .5))
rands <- data.frame(cbind(top = top, bot = bot, med = med, bin = seq(length(top))))
cohort <- ggplot(rands, aes(x = bin, y = med, ymin = bot, ymax = top))
cohort <- cohort + geom_pointrange()
cohort <- cohort + geom_hline(aes(yintercept = 0))
cohort <- cohort + labs(x = 'Time (My)', y = 'Cohort effect')
