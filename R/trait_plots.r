library(plyr)
library(reshape2)
library(ggplot2)
library(GGally)

source('../R/trait_analysis.r')

theme_set(theme_bw())
cbp <- c('#E69F00', '#56B4E9', '#009E73', '#F0E442', 
         '#0072B2', '#D55E00', '#CC79A7')

# total ~ oxygen
dfwin.bg <- lapply(win.bg, function(x) {
                   as.list(diff(rev(unlist(x))))})
tot.sum <- melt(dfwin.bg)
dfox <- diff(rev(oxyd))
tot.sum <- cbind(tot.sum, oxy = dfox[match(tot.sum$L2, names(dfox))])
tot.sum$L2 <- as.numeric(tot.sum$L2)
tot <- ggplot(tot.sum, aes(x = oxy, y = value))
tot <- tot + geom_point() + stat_smooth(method = 'lm')
tot <- tot + facet_wrap(~ L1, scales = 'free')
tot <- tot + labs(x = 'first differences oxygen',
                  y = 'first differences network statistic')
ggsave(file = '../doc/figure/tot_oxy.png', plot = tot)

# diet ~ oxygen
dfdt.bg <- lapply(dtwin.bg, function(x) {
                  lapply(x, function(y) {
                         as.list(diff(rev(unlist(y))))})})
dt.sum <- melt(dfdt.bg)
dt.sum <- cbind(dt.sum, oxy = dfox[match(dt.sum$L3, names(dfox))]) 
dt.sum$L3 <- as.numeric(dt.sum$L3)
dtox <- ggplot(dt.sum, aes(x = oxy, y = value, colour = L1))
dtox <- dtox + geom_point() + stat_smooth(method = 'lm', se = FALSE)
dtox <- dtox + facet_wrap(~ L2, scales = 'free')
dtox <- dtox + scale_color_manual(values = cbp)
dtox <- dtox + labs(x = 'first differences oxygen',
                    y = 'first differences network statistic')
ggsave(file = '../doc/figure/dt_oxy.png', plot = dtox)

# loco ~ oxygen
#dflf.bg <- lapply(lfwin.bg, function(x) {
#                  lapply(x, function(y) {
#                         as.list(diff(rev(unlist(y))))})})
#lf.sum <- melt(dflf.bg)
#lf.sum <- cbind(lf.sum, oxy = dfox[match(lf.sum$L3, names(dfox))]) 
#lf.sum$L3 <- as.numeric(lf.sum$L3)
#dtox <- ggplot(lf.sum, aes(x = oxy, y = value, colour = L1))
#dtox <- dtox + geom_point() + stat_smooth(method = 'lm', se = FALSE)
#dtox <- dtox + facet_wrap(~ L2, scales = 'free')
#dtox <- dtox + labs(x = 'first differences oxygen',
#                    y = 'first differences network statistic')
#ggsave(file = '../doc/figure/dt_oxy.png', plot = tot)
