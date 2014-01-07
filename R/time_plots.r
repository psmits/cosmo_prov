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

# just NA
nadt <- melt(dtwin.bg)
nadt$L3 <- as.numeric(nadt$L3)
ggnad <- ggplot(nadt, aes(x = L3, y = value, colour = L1))
ggnad <- ggnad + geom_line()
ggnad <- ggnad + stat_smooth(method = 'loess', se = FALSE, na.rm = TRUE)
ggnad <- ggnad + scale_color_manual(values = cbp,
                                    name = 'Dietary\nCategory')
ggnad <- ggnad + labs(x = 'Time (My)')
ggnad <- ggnad + facet_wrap(~ L2, nrow = 1, scales = 'free')
ggsave(file = '../doc/figure/na_dt.png', width = 15, height = 5, plot = ggnad)

# just Eur
erdt <- melt(dteur.bg)
erdt$L3 <- as.numeric(erdt$L3)
ggerd <- ggplot(erdt, aes(x = L3, y = value, colour = L1))
ggerd <- ggerd + geom_line()
ggerd <- ggerd + stat_smooth(method = 'loess', se = FALSE, na.rm = TRUE)
ggerd <- ggerd + scale_color_manual(values = cbp,
                                    name = 'Dietary\nCategory')
ggerd <- ggerd + labs(x = 'Time (My)')
ggerd <- ggerd + facet_wrap(~ L2, nrow = 1, scales = 'free')
ggsave(file = '../doc/figure/er_dt.png', width = 15, height = 5, plot = ggerd)


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

# just NA
nalf <- melt(lfwin.bg)
nalf$L3 <- as.numeric(nalf$L3)
ggnal <- ggplot(nalf, aes(x = L3, y = value, colour = L1))
ggnal <- ggnal + geom_line()
ggnal <- ggnal + stat_smooth(method = 'loess', se = FALSE, na.rm = TRUE)
ggnal <- ggnal + scale_color_manual(values = cbp,
                                    name = 'Locomotor\nCategory')
ggnal <- ggnal + labs(x = 'Time (My)')
ggnal <- ggnal + facet_wrap(~ L2, nrow = 1, scales = 'free')
ggsave(file = '../doc/figure/na_lf.png', width = 15, height = 5, plot = ggnal)

# just Eur
erlf <- melt(lfeur.bg)
erlf$L3 <- as.numeric(erlf$L3)
ggerl <- ggplot(erlf, aes(x = L3, y = value, colour = L1))
ggerl <- ggerl + geom_line()
ggerl <- ggerl + stat_smooth(method = 'loess', se = FALSE, na.rm = TRUE)
ggerl <- ggerl + scale_color_manual(values = cbp,
                                    name = 'Locomotor\nCategory')
ggerl <- ggerl + labs(x = 'Time (My)')
ggerl <- ggerl + facet_wrap(~ L2, nrow = 1, scales = 'free')
ggsave(file = '../doc/figure/er_lf.png', width = 15, height = 5, plot = ggerl)
