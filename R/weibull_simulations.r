library(rstan)
library(plyr)
library(MASS)
library(plyr)
library(parallel)
library(ggplot2)
library(grid)
library(reshape2)

theme_set(theme_bw())
cbp <- c('#E69F00', '#56B4E9', '#009E73', '#F0E442', 
         '#0072B2', '#D55E00', '#CC79A7')
theme_update(axis.text = element_text(size = 15),
             axis.title = element_text(size = 25),
             legend.text = element_text(size = 25),
             legend.title = element_text(size = 26),
             legend.key.size = unit(2, 'cm'),
             strip.text = element_text(size = 20))

set.seed(420)
n <- c(10, 100, 1000, 10000)
nn <- llply(n, function(x) rep(x, 1000))

sims <- llply(nn, function(x) {
              o <- llply(x, function(y) rweibull(y, shape = 1.3))
              o})

modes <- llply(sims, function(y) {
               aa <- llply(y, function(x) fitdistr(x, 'weibull'))
               aa <- laply(aa, function(x) x$estimate)
               aa})
names(modes) <- c('10', '100', '1000', '10000')
modes <- melt(modes)
modes$Var2 <- mapvalues(modes$Var2, 
                        from = unique(modes$Var2), 
                        to = c('alpha', 'sigma'))
vline <- data.frame(Var2 = c('alpha', 'sigma'), value = c(1.3, 1))


r <- ggplot(modes, aes(x = value))
r <- r + geom_vline(data = vline, aes(xintercept = value), colour = 'darkgrey')
r <- r + geom_histogram(aes(y = ..density..))
r <- r + facet_grid(Var2 ~ L1, labeller = label_parsed)
r <- r + labs(x = 'Estimated value', y = 'Prob. Density')
ggsave(r, file = '../doc/na_surv/figure/alpha_simulation.png', 
       width = 10, height = 8)
