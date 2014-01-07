library(mgcv)
library(reshape2)
library(ggplot2)
library(scales)

source('../R/na_mung.r')
source('../R/europe_mung.r')

source('../R/cosmo_prov.r')
source('../R/diet_dynamics.r')
source('../R/life_dynamics.r')
source('../R/oxygen_curve.r')

theme_set(theme_bw())
cbp <- c('#E69F00', '#56B4E9', '#009E73', '#F0E442', 
         '#0072B2', '#D55E00', '#CC79A7')

# zachos curve
ggzac <- ggplot(zachos, aes(x = age, y = -d18o)) 
ggzac <- ggzac + geom_point(aes(alpha = 0.3))
#ggzac <- ggzac + scale_color_manual(values = cbp)
ggzac <- ggzac + theme(legend.position = 'none')
ggzac <- ggzac + labs(x = 'Time (My)')
ggzac <- ggzac + stat_smooth(method = 'gam', 
                             formula = y ~ s(x, k = 13, bs = 'cs'))
ggzac <- ggzac + theme(axis.title = element_text(size = 19),
                       axis.text = element_text(size = 17))
ggsave(file = '../doc/figure/zachos.png', 
       width = 15, height = 10, plot = ggzac)

# standard bin
win.bg <- win.bg[1:4]
cont <- list(na = win.bg, eur = eurwin.bg)
bin.dat <- melt(cont)
bin.dat$L3 <- as.numeric(bin.dat$L3)
bin.dat <- bin.dat[!(bin.dat[, 1] == Inf | is.na(bin.dat[, 1])), ]
ggdat <- ggplot(bin.dat, aes(x = L3, y = value, linetype = L1)) + geom_line()
ggdat <- ggdat + labs(x = 'Time (My)')
ggdat <- ggdat + scale_linetype_manual(values = c(1,2),
                                       name = 'Region')
ggdat <- ggdat + stat_smooth(method = 'loess', se = FALSE, na.rm = TRUE)
ggdat <- ggdat + facet_wrap(~ L2, scales = 'free_y')
ggsave(file = '../doc/figure/gen_bin.png', 
       width = 15, height = 10, plot = ggdat)


# diet
dietdf <- list(na = dtwin.bg, eur = dteur.bg)
dietdf <- melt(dietdf)
dietdf$L4 <- as.numeric(dietdf$L4)
ggdiet <- ggplot(dietdf, aes(x = L4, y = value, colour = L2))
ggdiet <- ggdiet + geom_line()
ggdiet <- ggdiet + scale_color_manual(values = cbp)
ggdiet <- ggdiet + labs(x = 'Time (My)')
ggdiet <- ggdiet + facet_grid(L3 ~ L1, scales = 'free')
ggsave(file = '../doc/figure/diets.png', width = 10, height = 15, plot = ggdiet)


# locomotor
locodf <- list(na = lfwin.bg, eur = lfeur.bg)
locodf <- melt(locodf)
locodf$L4 <- as.numeric(locodf$L4)
ggloco <- ggplot(locodf, aes(x = L4, y = value, colour = L2))
ggloco <- ggloco + geom_line()
ggloco <- ggloco + scale_color_manual(values = cbp)
ggloco <- ggloco + labs(x = 'Time (My)')
ggloco <- ggloco + facet_grid(L3 ~ L1, scales = 'free')
ggsave(file = '../doc/figure/locos.png', width = 10, height = 15, plot = ggloco)
