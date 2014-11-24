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

# bioprovinces through time?

# posterior predictive checks
base.dur <- data.frame(dur = duration)
hist.sim <- melt(mm)  # simulations
hist.sim <- hist.sim[hist.sim$L1 %in% 1:12, ]
ppc.hist <- ggplot(hist.sim, aes(x = value))
ppc.hist <- ppc.hist + geom_histogram(data = base.dur, 
                                      aes(x = dur, y = ..density..),
                                      binwidth = 1, fill = 'grey', alpha = 0.5,
                                      colour = 'darkgrey')
ppc.hist <- ppc.hist + geom_histogram(aes(y = ..density..), binwidth = 1, 
                                      fill = 'blue', alpha = 0.4)
ppc.hist <- ppc.hist + labs(x = 'Duration', y = 'Prob. Density')
ppc.hist <- ppc.hist + facet_wrap( ~ L1, nrow = 3, ncol = 4)
#ggsave(ppc.hist, filename = '../doc/figure/histogram_ppc.png',
#       width = 15, height = 10)

means <- laply(mm, mean)
mean.dur <- mean(duration)
ppc.mean <- ggplot(data.frame(x = means), aes(x = x))
ppc.mean <- ppc.mean + geom_histogram(aes(y = ..density..), binwidth = 1)
ppc.mean <- ppc.mean + geom_vline(xintercept = mean.dur, 
                                  colour = 'blue', size = 2)
ppc.mean <- ppc.mean + labs(x = 'Duration', y = 'Prob. Density')
#ggsave(ppc.mean, filename = '../doc/figure/mean_ppc.png',
#       width = 15, height = 10)

# standardized residuals: mean and var? 
std.res <- melt(llply(mm, function(x) (duration - x) / sd(x)))
std.res <- std.res[std.res$L1 %in% 1:12, ]
std.res$index <- rep(seq(1:1921), 12)
ppc.res <- ggplot(std.res, aes(x = index, y = value))
ppc.res <- ppc.res + geom_point()
ppc.res <- ppc.res + facet_wrap( ~ L1, nrow = 3, ncol = 4)
#ggsave(ppc.res, filename = '../doc/figure/residual_plot.png',
#       width = 15, height = 10)


# survival function
emp.surv <- survfit(Surv(duration, extinct) ~ 1)
emp.surv <- data.frame(cbind(time = emp.surv$time, surv = emp.surv$surv))
emp.surv <- rbind(c(0, 1), emp.surv)

sim.surv <- llply(mm, function(x) survfit(Surv(x) ~ 1))
sim.surv <- llply(sim.surv, function(x) data.frame(cbind(time = x$time, 
                                                         surv = x$surv)))
sim.surv <- llply(sim.surv, function(x) rbind(c(0, 1), x))
sim.surv <- Reduce(rbind, Map(function(x, y) cbind(x, label = rep(y, nrow(x))),
                              x = sim.surv, y = seq(length(sim.surv))))
sim.surv <- sim.surv[sim.surv$label %in% 1:100, ]

soft <- ggplot(emp.surv, aes(x = time, y = surv))
soft <- soft + geom_step(data = sim.surv, aes(x = time, y = surv, group = label), 
                         colour = 'grey', alpha = 0.1)
soft <- soft + geom_step(size = 1, direction = 'hv')
soft <- soft + coord_cartesian(xlim = c(-0.5, max(duration) + 2))
#ggsave(soft, filename = '../doc/figure/survival_function.png',
#       width = 15, height = 10)

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
# better to compare as differences?
arb.eff <- base.inter$value
grd.eff <- base.inter$value + move.eff$value[1:ns]
scn.eff <- base.inter$value + move.eff$value[(ns + 1):(2 * ns)]
loceff <- melt(cbind(arb.eff, grd.eff, scn.eff))
loco <- ggplot(loceff, aes(x = value))
loco <- loco + geom_vline(xintercept = 0, colour = 'grey', size = 2)
loco <- loco + geom_histogram(aes(y = ..density..))
loco <- loco + facet_grid(Var2 ~ .)
loco <- loco + labs(x = 'Value', y = 'Prob. Density')
#ggsave(loco, filename = '../doc/figure/loco_est.png',
#       width = 15, height = 10)


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
diet <- diet + geom_vline(xintercept = 0, colour = 'grey', size = 2)
diet <- diet + geom_histogram(aes(y = ..density..))
diet <- diet + facet_grid(Var2 ~ .)
diet <- diet + labs(x = 'Value', y = 'Prob. Density')
#ggsave(diet, filename = '../doc/figure/diet_est.png',
#       width = 15, height = 10)

# effect of body size and occupancy
size.eff <- melted[melted$L1 == 'beta_size', 'value']  # base line
occ.eff <- melted[melted$L1 == 'beta_occ', 'value']  # base line
oth.eff <- melt(cbind(size.eff, occ.eff))
other <- ggplot(oth.eff, aes(x = value))
other <- other + geom_vline(xintercept = 0, colour = 'grey', size = 2)
other <- other + geom_histogram(aes(y = ..density..))
other <- other + facet_grid(Var2 ~ .)
other <- other + labs(x = 'Value', y = 'Prob. Density')
#ggsave(other, filename = '../doc/figure/other_est.png',
#       width = 15, height = 10)

# histogram of variance in cohort effect?
# for cohort effect, do point range with 80% quartile
top <- apply(mpost$rando, 2, function(x) quantile(x, .8))
bot <- apply(mpost$rando, 2, function(x) quantile(x, .2))
med <- apply(mpost$rando, 2, function(x) quantile(x, .5))
rands <- data.frame(cbind(top = top, bot = bot, med = med, 
                          bin = seq(length(top))))
rands$bin <- (rands$bin * 2) + 1
cohort <- ggplot(rands, aes(x = bin, y = med, ymin = bot, ymax = top))
cohort <- cohort + geom_hline(aes(yintercept = 0), colour = 'grey', size = 2)
cohort <- cohort + geom_pointrange()
cohort <- cohort + scale_x_continuous(breaks = seq(from = 0, to = 65, by = 5))
cohort <- cohort + labs(x = 'Time (My)', y = 'Cohort effect')
#ggsave(cohort, filename = '../doc/figure/cohort_est.png',
#       width = 15, height = 10)
