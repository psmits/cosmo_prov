library(ggplot2)
library(reshape2)
library(scales)
library(hexbin)
library(stringr)
library(grid)
library(survival)
library(GGally)
library(moments)
library(plyr)

nsim <- 100

pairwise.diffs <- function(x) {
  # create column combination pairs
  prs <- cbind(rep(1:ncol(x), each = ncol(x)), 1:ncol(x))
  col.diffs <- prs[prs[, 1] < prs[, 2], , drop = FALSE]
  # do pairwise differences 
  result <- x[, col.diffs[, 1]] - x[, col.diffs[, 2], drop = FALSE]
  # set colnames
  if(is.null(colnames(x)))
    colnames(x) <- 1:ncol(x)
  
  colnames(result) <- paste(colnames(x)[col.diffs[, 1]], ".vs.", 
                            colnames(x)[col.diffs[, 2]], sep = "")
  result
}



theme_set(theme_bw())
cbp <- c('#E69F00', '#56B4E9', '#009E73', '#F0E442', 
         '#0072B2', '#D55E00', '#CC79A7')
theme_update(axis.text = element_text(size = 20),
             axis.title = element_text(size = 30),
             legend.text = element_text(size = 25),
             legend.title = element_text(size = 26),
             legend.key.size = unit(2, 'cm'),
             strip.text = element_text(size = 20))


# posterior predictive checks
base.dur <- data.frame(dur = duration)
hist.sim <- melt(pm)  # simulations
hist.sim <- hist.sim[hist.sim$L1 %in% 1:12, ]
ppc.hist <- ggplot(hist.sim, aes(x = value))
ppc.hist <- ppc.hist + geom_histogram(data = base.dur, 
                                      aes(x = dur, y = ..density..),
                                      fill = 'grey', alpha = 0.5,
                                      colour = 'darkgrey')
ppc.hist <- ppc.hist + geom_histogram(aes(y = ..density..), binwidth = .2, 
                                      fill = 'blue', alpha = 0.4)
ppc.hist <- ppc.hist + labs(x = 'Duration', y = 'Prob. Density')
ppc.hist <- ppc.hist + facet_wrap( ~ L1, nrow = 3, ncol = 4)
ggsave(ppc.hist, filename = '../doc/na_surv/figure/histogram_ppc.png',
       width = 15, height = 10)

means <- laply(pm, mean)
mean.dur <- mean(duration)
ppc.mean <- ggplot(data.frame(x = means), aes(x = x))
ppc.mean <- ppc.mean + geom_histogram(aes(y = ..density..), binwidth = .2)
ppc.mean <- ppc.mean + geom_vline(xintercept = mean.dur, 
                                  colour = 'blue', size = 2)
ppc.mean <- ppc.mean + labs(x = 'Duration', y = 'Prob. Density')
ggsave(ppc.mean, filename = '../doc/na_surv/figure/mean_ppc.png',
       width = 10, height = 10)

quant <- laply(pm, function(x) c(mean = mean(x), quantile(x, c(.25, .5, .75))))
quant <- melt(quant)
quant.dur <- c(mean = mean(duration), quantile(duration, c(.25, .5, .75)))
quant.dur <- melt(quant.dur)
quant.dur$Var2 <- rownames(quant.dur)

ppc.quant <- ggplot(quant, aes(x = value))
ppc.quant <- ppc.quant + geom_histogram(aes(y = ..density..), binwidth = .2)
ppc.quant <- ppc.quant + geom_vline(data = quant.dur, aes(xintercept = value), 
                                    colour = 'blue', size = 2)
ppc.quant <- ppc.quant + labs(x = 'Duration', y = 'Prob. Density')
ppc.quant <- ppc.quant + facet_wrap(~ Var2, ncol = 2)
ggsave(ppc.quant, filename = '../doc/na_surv/figure/quant_ppc.png',
       width = 10, height = 10)

# deviance residuals
std.res <- melt(pm.res)
std.res <- std.res[std.res$L1 %in% 1:12, ]
std.res$index <- rep(seq(1:1921), 12)
ppc.res <- ggplot(std.res, aes(x = index, y = value))
ppc.res <- ppc.res + geom_hline(aes(yintercept = 0), 
                                colour = 'darkgrey', size = 1)
ppc.res <- ppc.res + geom_hline(aes(yintercept = 2), 
                                colour = 'darkgrey', size = 1, 
                                linetype = 'dashed')
ppc.res <- ppc.res + geom_hline(aes(yintercept = -2), 
                                colour = 'darkgrey', size = 1, 
                                linetype = 'dashed')
ppc.res <- ppc.res + geom_point(alpha = 0.5, size = 1)
ppc.res <- ppc.res + facet_wrap( ~ L1, nrow = 3, ncol = 4)
ppc.res <- ppc.res + labs(x = '', y = 'deviance residuals')
ppc.res <- ppc.res + theme(axis.ticks.x = element_blank(),
                           axis.text.x = element_blank())
ggsave(ppc.res, filename = '../doc/na_surv/figure/residual_plot.png',
       width = 15, height = 10)

skew.res <- laply(pm.res, moments::skewness)
var.res <- laply(pm.res, function(x) var(x))
res.sum <- melt(cbind(skew.res, var.res))
res.sum$Var2 <- as.character(res.sum$Var2)
res.sum$Var2[res.sum$Var2 == 'skew.res'] <- 'residual skewness'
res.sum$Var2[res.sum$Var2 == 'var.res'] <- 'residual variance'
ppc.sum <- ggplot(res.sum, aes(x = value))
ppc.sum <- ppc.sum + geom_vline(xintercept = 0, colour = 'grey', size = 2)
ppc.sum <- ppc.sum + geom_histogram(aes(y = ..density..), binwidth = 0.2)
ppc.sum <- ppc.sum + facet_grid(Var2 ~ .)
ppc.sum <- ppc.sum + labs(x = 'value', y = 'Prob. Density')
ggsave(ppc.sum, filename = '../doc/na_surv/figure/res_sum_plot.png',
       width = 5, height = 10)


# survival function
condition <- extinct
condition[extinct == 1 & duration == 1] <- 2
emp.surv <- survfit(Surv(time = duration, time2 = duration, 
                         event = condition, type = 'interval') ~ 1)
emp.surv <- data.frame(cbind(time = emp.surv$time, surv = emp.surv$surv))
emp.surv <- rbind(c(0, 1), emp.surv)

sim.surv <- llply(pm, function(x) survfit(Surv(x) ~ 1))
sim.surv <- llply(sim.surv, function(x) data.frame(cbind(time = x$time, 
                                                         surv = x$surv)))
sim.surv <- llply(sim.surv, function(x) rbind(c(0, 1), x))
sim.surv <- Reduce(rbind, Map(function(x, y) cbind(x, label = rep(y, nrow(x))),
                              x = sim.surv, y = seq(length(sim.surv))))
sim.surv <- sim.surv[sim.surv$label %in% 1:nsim, ]

exp.surv <- llply(ee, function(x) survfit(Surv(x) ~ 1))
exp.surv <- llply(exp.surv, function(x) data.frame(cbind(time = x$time, 
                                                         surv = x$surv)))
exp.surv <- llply(exp.surv, function(x) rbind(c(0, 1), x))
exp.surv <- Reduce(rbind, Map(function(x, y) cbind(x, label = rep(y, nrow(x))),
                              x = exp.surv, y = seq(length(exp.surv))))
exp.surv <- exp.surv[exp.surv$label %in% 1:nsim, ]

mix.surv <- rbind(cbind(sim.surv, lab = rep('Weibull', nrow(sim.surv))),
                  cbind(exp.surv, lab = rep('Exponential', nrow(exp.surv))))

soft <- ggplot(emp.surv, aes(x = time, y = surv))
soft <- soft + geom_step(data = mix.surv, aes(x = time, y = surv, group = label), 
                         colour = 'grey', alpha = .2)
soft <- soft + geom_step(size = 1, direction = 'hv')
soft <- soft + coord_cartesian(xlim = c(-0.5, max(duration) + 2))
soft <- soft + facet_grid(. ~ lab)
soft <- soft + labs(x = 'Duration', y = 'P(T > t)')
ggsave(soft, filename = '../doc/na_surv/figure/survival_function.png',
       width = 15, height = 10)



# marginal posteriors
melted <- melt(mpost)
scale.melted <- melt(phypost)
ns <- length(phypost$lp__)

# create the different combinations
regcoef <- scale.melted[str_detect(scale.melted$L1, 'beta'), ]
base.inter <- scale.melted[scale.melted$L1 == 'beta_inter', ]  # base line
move.eff <- scale.melted[scale.melted$L1 == 'beta_move', ]  # move effect
diet.eff <- scale.melted[scale.melted$L1 == 'beta_diet', ]  # diet effect

# make fancy versions of these
# lower values, lower risk
# these are particularilty useful for summary statistics
# effect of locomotor category
arb <- base.inter$value
grd <- base.inter$value + move.eff$value[1:ns]
scn <- base.inter$value + move.eff$value[(ns + 1):(2 * ns)]
loceff <- melt(cbind(arb, grd, scn))
loco <- ggplot(loceff, aes(x = value))
loco <- loco + geom_vline(xintercept = 0, colour = 'grey', size = 2)
loco <- loco + geom_histogram(aes(y = ..density..))
loco <- loco + facet_grid(Var2 ~ .)
loco <- loco + labs(x = 'Value', y = 'Prob. Density')
ggsave(loco, filename = '../doc/na_surv/figure/loco_est.png',
       width = 10, height = 10)

# better to compare as differences?
loco.diff <- melt(pairwise.diffs(cbind(arb, grd, scn)))
lodf <- ggplot(loco.diff, aes(x = value))
lodf <- lodf + geom_vline(xintercept = 0, colour = 'grey', size = 2)
lodf <- lodf + geom_histogram(aes(y = ..density..))
lodf <- lodf + facet_grid(Var2 ~ .)
lodf <- lodf + labs(x = 'Value', y = 'Prob. Density')
ggsave(lodf, filename = '../doc/na_surv/figure/loco_diff_est.png',
       width = 10, height = 10)


# effect of dietary category
crn <- base.inter$value
hrb <- base.inter$value + diet.eff$value[1:ns]
ist <- base.inter$value + diet.eff$value[(ns + 1):(2 * ns)]
omn <- base.inter$value + diet.eff$value[(2 * ns + 1):(3 * ns)]
deteff <- melt(cbind(crn, hrb, ist, omn))
diet <- ggplot(deteff, aes(x = value))
diet <- diet + geom_vline(xintercept = 0, colour = 'grey', size = 2)
diet <- diet + geom_histogram(aes(y = ..density..))
diet <- diet + facet_grid(Var2 ~ .)
diet <- diet + labs(x = 'Value', y = 'Prob. Density')
ggsave(diet, filename = '../doc/na_surv/figure/diet_est.png',
       width = 10, height = 10)

# better to compare as differences?
diet.diff <- melt(pairwise.diffs(cbind(crn, hrb, ist, omn)))
didf <- ggplot(diet.diff, aes(x = value))
didf <- didf + geom_vline(xintercept = 0, colour = 'grey', size = 2)
didf <- didf + geom_histogram(aes(y = ..density..))
didf <- didf + facet_grid(Var2 ~ .)
didf <- didf + labs(x = 'Value', y = 'Prob. Density')
ggsave(didf, filename = '../doc/na_surv/figure/diet_diff_est.png',
       width = 10, height = 10)


# effect of body size and occupancy
size.eff <- melted[melted$L1 == 'beta_size', 'value']  # base line
occ.eff <- melted[melted$L1 == 'beta_occ', 'value']  # base line
oth.eff <- melt(cbind(size.eff = size.eff[1:2000], 
                      occ.eff = occ.eff[1:2000]))
oth.eff$label <- rep('unstandardized', nrow(oth.eff))

sc.size.eff <- scale.melted[scale.melted$L1 == 'beta_size', 'value']
sc.occ.eff <- scale.melted[scale.melted$L1 == 'beta_occ', 'value']
sc.oth.eff <- melt(cbind(sc.size.eff, sc.occ.eff))
sc.oth.eff$label <- rep('standardized', nrow(sc.oth.eff))

sc.oth.eff$Var2 <- oth.eff$Var2
oth.eff <- rbind(oth.eff, sc.oth.eff)

other <- ggplot(oth.eff, aes(x = value))
other <- other + geom_vline(xintercept = 0, colour = 'grey', size = 2)
other <- other + geom_histogram(aes(y = ..density..))
other <- other + facet_grid(Var2 ~ label)
other <- other + labs(x = 'Value', y = 'Prob. Density')
ggsave(other, filename = '../doc/na_surv/figure/other_est.png',
       width = 10, height = 10)


# variance partition coefficient
s.y <- var.star[, 1]
s.c <- var.star[, 3]
s.p <- var.star[, 2]

indiv.part <- s.y / (s.y + s.c + s.p)
cohort.part <- s.c / (s.y + s.c + s.p)
phylo.part <- s.p / (s.y + s.c + s.p)

var.parts <- melt(cbind(individual = indiv.part, 
                        cohort = cohort.part, 
                        phylogeny = phylo.part))
gvar <- ggplot(var.parts, aes(x = value))
gvar <- gvar + geom_histogram(aes(y = ..density..))
gvar <- gvar + facet_grid(Var2 ~ .)
gvar <- gvar + scale_x_continuous(limits = c(0, 1))
gvar <- gvar + labs(x = 'Variance partition\ncoefficient', y = 'Prob. Density')
ggsave(gvar, filename = '../doc/na_surv/figure/variance_est.png',
       width = 5, height = 10)


# for cohort effect, do point range with 80% quartile
top <- apply(phypost$rando, 2, function(x) quantile(x, .8))
bot <- apply(phypost$rando, 2, function(x) quantile(x, .2))
med <- apply(phypost$rando, 2, function(x) quantile(x, .5))
rands <- data.frame(cbind(top = top, bot = bot, med = med, 
                          bin = seq(length(top))))
rands$bin <- (rands$bin * 2) + 1
cohort <- ggplot(rands, aes(x = bin, y = med, ymin = bot, ymax = top))
cohort <- cohort + geom_hline(aes(yintercept = 0), colour = 'grey', size = 2)
cohort <- cohort + geom_pointrange()
cohort <- cohort + scale_x_reverse(breaks = seq(from = 0, to = 65, by = 5))
cohort <- cohort + labs(x = 'Time (My)', y = 'Cohort effect')
ggsave(cohort, filename = '../doc/na_surv/figure/cohort_est.png',
       width = 15, height = 10)


# estimate of alpha
wei.haz <- function(time, scale, alpha) {
  # lambda is the inverse-scale parameter
  tt <- time**(alpha - 1)
  scale * alpha * tt
}

xx <- seq(0, 12, by = 0.01)
wz <- list()
for(ii in seq(nsim)) {
  a <- sample(phypost$alpha, 1)
  wz[[ii]] <- wei.haz(xx, exp(-(sample(phypost$beta_inter, 1)) / a ), a)
}

ma <- median(phypost$alpha)
mid.haz <- wei.haz(xx, exp(-(median(phypost$beta_inter) / ma)), ma)
mid.haz <- data.frame(time = xx, hazard = mid.haz)

wzm <- llply(wz, function(x) data.frame(time = xx, hazard = x))
wzm <- Reduce(rbind, Map(function(x, y) cbind(x, label = rep(y, nrow(x))), 
                         x = wzm, y = seq(length(wzm))))
haz <- ggplot(mid.haz, aes(x = time, y = hazard))
haz <- haz + geom_line(data = wzm, aes(x = time, y = hazard, group = label),
                       colour = 'darkgrey', alpha = 0.2)
haz <- haz + geom_line(colour = 'blue')
haz <- haz + labs(x = 'Duration', y = 'h(t)')
ggsave(haz, filename = '../doc/na_surv/figure/haz_est.png',
       width = 15, height = 10)
