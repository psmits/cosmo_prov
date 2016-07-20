library(ggplot2)
library(gridExtra)
library(reshape2)
library(scales)
library(hexbin)
library(stringr)
library(grid)
library(survival)
library(moments)
library(plyr)
source('../R/waic.r')
source('../R/multiplot.r')
source('../R/surv_post_sim.r')
wei.waic <- waic(phy.scalemfit)
#exp.waic <- waic(escalefit)

nsim <- 1000
duration <- c(data$dur_unc, data$dur_cen)

pairwise.diffs <- function(x) {
  # create column combination pairs
  prs <- cbind(rep(1:ncol(x), each = ncol(x)), 1:ncol(x))
  col.diffs <- prs[prs[, 1] < prs[, 2], , drop = FALSE]
  # do pairwise differences 
  result <- x[, col.diffs[, 1]] - x[, col.diffs[, 2], drop = FALSE]
  # set colnames
  if(is.null(colnames(x)))
    colnames(x) <- 1:ncol(x)

  colnames(result) <- paste('beta[', colnames(x)[col.diffs[, 1]], "] - beta[", 
                            colnames(x)[col.diffs[, 2]], ']', sep = "")
  result
}


theme_set(theme_bw())
cbp <- c('#E69F00', '#56B4E9', '#009E73', '#F0E442', 
         '#0072B2', '#D55E00', '#CC79A7')
theme_update(axis.text = element_text(size = 6),
             axis.title = element_text(size = 10),
             legend.text = element_text(size = 8),
             legend.title = element_text(size = 9),
             legend.key.size = unit(0.75, 'cm'),
             strip.text = element_text(size = 8))


# posterior predictive checks
quant <- laply(pm, function(x) c(mean = mean(x), quantile(x, c(.25, .5, .75))))
quant <- melt(quant)
quant.dur <- c(mean = mean(duration), quantile(duration, c(.25, .5, .75)))
quant.dur <- melt(quant.dur)
quant.dur$Var2 <- rownames(quant.dur)

# all four of the major point checks
ppc.quant <- ggplot(quant, aes(x = value))
ppc.quant <- ppc.quant + geom_histogram(aes(y = ..density..), binwidth = .2)
ppc.quant <- ppc.quant + geom_vline(data = quant.dur, aes(xintercept = value), 
                                    colour = 'blue', size = 1)
ppc.quant <- ppc.quant + labs(x = 'Duration (2 My bins)', y = 'Prob. Density')
ppc.quant <- ppc.quant + facet_wrap(~ Var2, ncol = 2)
ggsave(ppc.quant, filename = '../doc/na_surv/figure/quant_ppc_noj.pdf',
       width = 3.42, height = 2.75, dpi = 750)
#ggsave(ppc.quant, filename = '../doc/na_surv/figure/quant_ppc.pdf',
#       width = 3.42, height = 2.75, dpi = 750)
ggsave(ppc.quant, filename = '../doc/na_surv/figure/quant_ppc_pres_noj.pdf',
       width = 4.7, height = 3.5, dpi = 750)
#ggsave(ppc.quant, filename = '../doc/na_surv/figure/quant_ppc_pres.pdf',
#       width = 4.7, height = 3.5, dpi = 750)

# deviance residuals
# change this to be x = duration, y = residual
std.res <- melt(pm.res)
std.res <- std.res[std.res$L1 %in% 1:12, ]
std.res$index <- rep(duration, 12)
ppc.res <- ggplot(std.res, aes(x = index, y = value))
ppc.res <- ppc.res + geom_hline(aes(yintercept = 0), 
                                colour = 'darkgrey', size = 1)
ppc.res <- ppc.res + geom_hline(aes(yintercept = 2), 
                                colour = 'darkgrey', size = 1, 
                                linetype = 'dashed')
ppc.res <- ppc.res + geom_hline(aes(yintercept = -2), 
                                colour = 'darkgrey', size = 1, 
                                linetype = 'dashed')
ppc.res <- ppc.res + geom_point(alpha = 0.1, size = 0.5, position = 'jitter')
ppc.res <- ppc.res + facet_wrap( ~ L1, nrow = 3, ncol = 4)
ppc.res <- ppc.res + labs(x = 'Duration (2 My bins)', y = 'Deviance residual')
ggsave(ppc.res, filename = '../doc/na_surv/figure/residual_plot_noj.pdf',
       width = 3.42, height = 2.25, dpi = 750)
#ggsave(ppc.res, filename = '../doc/na_surv/figure/residual_plot.pdf',
#       width = 3.42, height = 2.25, dpi = 750)


# survival function
condition <- nojanis.extinct
#condition <- extinct
condition[nojanis.extinct == 1 & nojanis.duration == 1] <- 2
#condition[extinct == 1 & duration == 1] <- 2
emp.surv <- survfit(Surv(time = duration, time2 = duration, 
                         event = condition, type = 'interval') ~ 1)
emp.surv <- data.frame(cbind(time = emp.surv$time, surv = emp.surv$surv))
emp.surv <- rbind(c(0, 1), emp.surv)

sim.surv <- llply(pm, function(x) survfit(Surv(x) ~ 1))
sim.surv <- llply(sim.surv, function(x) {
                  y <- data.frame(cbind(time = x$time, surv = x$surv))
                  y <- rbind(c(0, 1), y)
                  y})
sim.surv <- Reduce(rbind, Map(function(x, y) {
                              x$group <- y
                              x},
                              x = sim.surv, y = seq(length(sim.surv))))
sim.surv <- sim.surv[sim.surv$group %in% 1:nsim, ]
#sim.surv$lab <- 'Weibull'

#exp.surv <- llply(ee, function(x) survfit(Surv(x) ~ 1))
#exp.surv <- llply(exp.surv, function(x) {
#                  y <- data.frame(cbind(time = x$time, surv = x$surv))
#                  y <- rbind(c(0, 1), y)
#                  y})
#exp.surv <- llply(exp.surv, function(x) rbind(c(0, 1), x))
#exp.surv <- Reduce(rbind, Map(function(x, y) {
#                              x$group <- y
#                              x},
#                              x = exp.surv, y = seq(length(exp.surv))))
#exp.surv <- exp.surv[exp.surv$group %in% 1:nsim, ]
#exp.surv$lab <- 'Exponential'
mix.surv <- sim.surv

soft <- ggplot(emp.surv, aes(x = time, y = surv))
soft <- soft + geom_line(data = mix.surv, aes(x = time, y = surv, 
                                              group = group),
                         colour = 'grey', alpha = .2)
soft <- soft + geom_line(size = 1)
soft <- soft + coord_cartesian(xlim = c(-0.5, max(duration) + 2))
#soft <- soft + facet_grid(. ~ lab)
soft <- soft + labs(x = 'Duration (2 My bins)', y = 'P(T > t)')
ggsave(soft, filename = '../doc/na_surv/figure/survival_function_noj.tiff',
       width = 3.42, height = 3.42, dpi = 750)
ggsave(soft, filename = '../doc/na_surv/figure/survival_function_pres_noj.pdf',
       width = 4.7, height = 3.5, dpi = 750)
#ggsave(soft, filename = '../doc/na_surv/figure/survival_function.tiff',
#       width = 3.42, height = 3.42, dpi = 750)
#ggsave(soft, filename = '../doc/na_surv/figure/survival_function_pres.pdf',
#       width = 4.7, height = 3.5, dpi = 750)

# do each cohort


# marginal posteriors
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

# better to compare as differences?
loco.diff <- melt(pairwise.diffs(cbind(arb, grd, scn)))
relab.x <- scale_x_discrete(labels = 
                            c('beta[arb] - beta[grd]' = 
                              expression(beta[arb]-beta[grd]), 
                              'beta[arb] - beta[scn]' = 
                              expression(beta[arb]-beta[scn]), 
                              'beta[grd] - beta[scn]' = 
                              expression(beta[grd]-beta[scn])))
lodf <- ggplot(loco.diff, aes(x = Var2, y = value))
lodf <- lodf + geom_hline(yintercept = 0, colour = 'grey', size = 1)
lodf <- lodf + geom_violin()# + geom_boxplot(width = 0.1, outlier.size = 1)
lodf <- lodf + relab.x + theme(axis.text.x = element_text(size = 7))
lodf <- lodf + labs(x = '', y = 'Estimated difference', title = 'A')
lodf <- lodf + theme(plot.title = element_text(hjust = 0, size = 10))
ggsave(lodf, filename = '../doc/na_surv/figure/loco_diff_est_noj.pdf',
       width = 3.42, height = 2.75, dpi = 750)
#ggsave(lodf, filename = '../doc/na_surv/figure/loco_diff_est.pdf',
#       width = 3.42, height = 2.75, dpi = 750)


# effect of dietary category
crn <- base.inter$value
hrb <- base.inter$value + diet.eff$value[1:ns]
ist <- base.inter$value + diet.eff$value[(ns + 1):(2 * ns)]
omn <- base.inter$value + diet.eff$value[(2 * ns + 1):(3 * ns)]

# better to compare as differences?
diet.diff <- melt(pairwise.diffs(cbind(crn, hrb, ist, omn)))
relab.x <- scale_x_discrete(labels = 
                            c('beta[crn] - beta[hrb]' = 
                              expression(beta[crn]-beta[hrb]), 
                              'beta[crn] - beta[ist]' = 
                              expression(beta[crn]-beta[ist]), 
                              'beta[crn] - beta[omn]' =
                              expression(beta[crn]-beta[omn]),
                              'beta[hrb] - beta[ist]' =
                              expression(beta[hrb]-beta[ist]),
                              'beta[hrb] - beta[omn]' = 
                              expression(beta[hrb]-beta[omn]),
                              'beta[ist] - beta[omn]' =
                              expression(beta[ist]-beta[omn])))
didf <- ggplot(diet.diff, aes(x = Var2, y = value))
didf <- didf + geom_hline(yintercept = 0, colour = 'grey', size = 1)
didf <- didf + geom_violin()# + geom_boxplot(width = 0.1, outlier.size = 1)
didf <- didf + relab.x + theme(axis.text.x = element_text(size = 7))
didf <- didf + labs(x = '', y = 'Estimated difference', title = 'B')
didf <- didf + theme(plot.title = element_text(hjust = 0, size = 10))
ggsave(didf, filename = '../doc/na_surv/figure/diet_diff_est_noj.pdf',
       width = 3.42, height = 2.75, dpi = 750)
#ggsave(didf, filename = '../doc/na_surv/figure/diet_diff_est.pdf',
#       width = 3.42, height = 2.75, dpi = 750)

tiff('../doc/na_surv/figure/trait_eff_noj.tiff', 
     width = 3.42, height = 5.5, units = 'in', res = 750)
multiplot(lodf, didf, cols = 1)
dev.off()
#tiff('../doc/na_surv/figure/trait_eff.tiff', 
#     width = 3.42, height = 5.5, units = 'in', res = 750)
#multiplot(lodf, didf, cols = 1)
#dev.off()


# effect of body size and occupancy
sc.size.eff <- scale.melted[scale.melted$L1 == 'beta_size', 'value']
sc.occ.eff <- scale.melted[scale.melted$L1 == 'beta_occ', 'value']
sc.oth.eff <- melt(cbind(sc.size.eff, sc.occ.eff))
sc.oth.eff$Var2 <- as.character(sc.oth.eff$Var2)
sc.oth.eff$Var2 <- mapvalues(sc.oth.eff$Var2, c('sc.occ.eff', 'sc.size.eff'), 
                             c('beta[occupancy]', 'beta[size]'))

sum(sc.size.eff < 0) / length(sc.size.eff)
mean(sc.size.eff)
sd(sc.size.eff)
sum(sc.occ.eff < 0) / length(sc.occ.eff)
mean(sc.occ.eff)
sd(sc.occ.eff)

other <- ggplot(sc.oth.eff, aes(x = value))
other <- other + geom_vline(xintercept = 0, colour = 'grey', size = 1)
other <- other + geom_histogram(aes(y = ..density..))
other <- other + facet_grid(Var2 ~ . , labeller = label_parsed)
other <- other + labs(x = 'Parameter estimate', y = 'Prob. Density')
ggsave(other, filename = '../doc/na_surv/figure/other_est_noj.pdf',
       width = 3.42, height = 2.00, dpi = 750)
#ggsave(other, filename = '../doc/na_surv/figure/other_est.pdf',
#       width = 3.42, height = 2.00, dpi = 750)


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
gvar <- gvar + scale_x_continuous(limits = c(0, 1))
gvar.v <- gvar + facet_grid(Var2 ~ .)
gvar.v <- gvar.v + labs(x = 'Variance partition\ncoefficient', 
                        y = 'Prob. Density')
ggsave(gvar.v, filename = '../doc/na_surv/figure/variance_est_noj.pdf',
       width = 3.42, height = 3.00, dpi = 750)
#ggsave(gvar.v, filename = '../doc/na_surv/figure/variance_est.pdf',
#       width = 3.42, height = 3.00, dpi = 750)
gvar.h <- gvar + coord_flip() + facet_grid(~ Var2)
gvar.h <- gvar.h + labs(x = 'Variance partition coefficient', 
                        y = 'Prob. Density')
ggsave(gvar.h, filename = '../doc/na_surv/figure/variance_est_pres_noj.pdf',
       width = 4.7, height = 3.5, dpi = 750)
#ggsave(gvar.h, filename = '../doc/na_surv/figure/variance_est_pres.pdf',
#       width = 4.7, height = 3.5, dpi = 750)


# for cohort effect, do point range with 80% quartile
top <- apply(phypost$rando, 2, function(x) quantile(x, .8))
bot <- apply(phypost$rando, 2, function(x) quantile(x, .2))
med <- apply(phypost$rando, 2, function(x) quantile(x, .5))
rands <- data.frame(cbind(top = top, bot = bot, med = med, 
                          bin = seq(length(top))))
rands$bin <- (rands$bin * 2) + 1
cohort <- ggplot(rands, aes(x = bin, y = med, ymin = bot, ymax = top))
cohort <- cohort + geom_hline(aes(yintercept = 0), colour = 'grey', size = 1)
cohort <- cohort + geom_pointrange(size = 0.75)
cohort <- cohort + scale_x_reverse(breaks = seq(from = 0, to = 65, by = 5))
cohort <- cohort + labs(x = 'Time (Mya)', y = 'Estimated cohort effect')
ggsave(cohort, filename = '../doc/na_surv/figure/cohort_est_noj.pdf',
       width = 3.42, height = 2.25, dpi = 750)
ggsave(cohort, filename = '../doc/na_surv/figure/cohort_est_pres_noj.pdf',
       width = 4.7, height = 3.5, dpi = 750)
#ggsave(cohort, filename = '../doc/na_surv/figure/cohort_est.pdf',
#       width = 3.42, height = 2.25, dpi = 750)
#ggsave(cohort, filename = '../doc/na_surv/figure/cohort_est_pres.pdf',
#       width = 4.7, height = 3.5, dpi = 750)


# estimate of alpha
wei.haz <- function(time, scale, alpha) {
  # lambda is the inverse-scale parameter
  tt <- (time / scale)**(alpha - 1)
  (alpha / scale) * tt
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
haz <- haz + labs(x = 'Duration (2 My bins)', y = 'h(t)')
ggsave(haz, filename = '../doc/na_surv/figure/haz_est_noj.pdf',
       width = 3.42, height = 2.25, dpi = 750)
#ggsave(haz, filename = '../doc/na_surv/figure/haz_est.pdf',
#       width = 3.42, height = 2.25, dpi = 750)
