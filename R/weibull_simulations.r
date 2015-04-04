library(rstan)
library(plyr)
library(MASS)
library(plyr)
library(parallel)
library(ggplot2)
library(reshape2)

mod <- stan(file = '../stan/shape_est.stan')
set.seed(420)
n <- c(10, 100, 1000, 10000)
nn <- llply(n, function(x) rep(x, 100))

sims <- llply(nn, function(x) {
              o <- llply(x, function(y) rweibull(y, shape = 1.3))
              o})

fits <- list()
for(ii in seq(length(sims))) {
  ins <- list()
  for(jj in seq(length(sims[[ii]]))) {
    fit.list <- mclapply(1:4, mc.cores = detectCores(),
                         function(x) stan(fit = mod,
                                          data = list(N = n[ii], 
                                                      y = sims[[ii]][[jj]]),
                                          seed = 420,
                                          chains = 1, chain_id = x,
                                          refresh = -1))
    ins[[jj]] <- sflist2stanfit(fit.list)
  }
  fits[[ii]] <- ins
}

means <- llply(fits, function(y) 
               Reduce(rbind, llply(y, function(x) 
                                   laply(extract(x, permuted = TRUE), mean))))
means <- llply(means, function(x) {
               x <- data.frame(x)
               names(x) <- c('alpha', 'sigma', 'lp__')
               x[, 1:2]})
means <- llply(means, function(x) melt(x))
means <- Reduce(rbind, Map(function(x, y) 
                           data.frame(x, size = rep(paste0('n == ', y), 
                                                    length(x))), 
                           x = means, y = n))
vline <- data.frame(variable = c('alpha', 'sigma'), value = c(1.3, 1))
r <- ggplot(means, aes(x = value))
r <- r + geom_vline(data = vline, aes(xintercept = value), colour = 'darkgrey')
r <- r + geom_histogram(aes(y = ..density..))
r <- r + facet_grid(variable ~ size, labeller = label_parsed)
ggsave(r, file = '../doc/na_surv/figure/alpha_simulation.png', 
       width = 3.5, height = 4)
