library(mgcv)
library(reshape2)
library(ggplot2)
library(scales)

source('../R/cosmo_prov.r')
#source('../R/diet_dynamics.r')
#source('../R/life_dynamics.r')
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
cont <- list(na = win.bg, eur = eurwin.bg)
bin.dat <- melt(cont)
bin.dat$L3 <- as.numeric(bin.dat$L3)
bin.dat <- bin.dat[!(bin.dat[, 1] == Inf | is.na(bin.dat[, 1])), ]
ggdat <- ggplot(bin.dat, aes(x = L3, y = value, colour = L1)) + geom_line()
ggdat <- ggdat + facet_wrap(~ L2, scales = 'free_y')
ggdat <- ggdat + scale_color_manual(values = cbp)
#ggdat <- ggdat + theme(legend.position = 'none')
ggdat <- ggdat + labs(x = 'Time (My)')
ggdat <- ggdat + stat_smooth(method = 'loess', se = FALSE, na.rm = TRUE)
ggsave(file = '../doc/figure/gen_bin.png', 
       width = 15, height = 10, plot = ggdat)


# diet
# na diet
bin.diet <- melt(dtwin.bg)
bin.diet$L3 <- as.numeric(bin.diet$L3)
nadiet <- ggplot(bin.diet, aes(x = L3, y = value, colour = L1))
nadiet <- nadiet + geom_line(alpha = 0.5)
nadiet <- nadiet + facet_wrap(~ L2, scales = 'free_y')
nadiet <- nadiet + scale_color_manual(values = cbp)
nadiet <- nadiet + labs(x = 'Time (My)')
nadiet <- nadiet + theme(legend.position = 'none', 
                         axis.title = element_text(size = 19),
                         axis.text = element_text(size = 17))
nadiet <- nadiet + stat_smooth(method = 'loess', se = FALSE)
ggsave(file = '../doc/figure/diet_bin.png', 
       width = 15, height = 10, plot = nadiet)

# europe diet
eub.diet <- melt(dteur.bg)
eub.diet$L3 <- as.numeric(eub.diet$L3)
eub.diet <- eub.diet[!(eub.diet[, 1] == Inf | is.na(eub.diet[, 1])), ]
eudiet <- ggplot(eub.diet, aes(x = L3, y = value, colour = L1))
eudiet <- eudiet + geom_line(alpha = 0.5)
eudiet <- eudiet + facet_wrap(~ L2, scales = 'free_y')
eudiet <- eudiet + scale_color_manual(values = cbp)
eudiet <- eudiet + labs(x = 'Time (My)')
eudiet <- eudiet + theme(legend.position = 'none', 
                         axis.title = element_text(size = 19),
                         axis.text = element_text(size = 17))
eudiet <- eudiet + stat_smooth(method = 'loess', se = FALSE)
ggsave(file = '../doc/figure/eudt_bin.png', 
       width = 15, height = 10, plot = nadiet)

# locomotor
bin.loco <- melt(lfwin.bg)
bin.loco$L3 <- as.numeric(bin.loco$L3)
ggloco <- ggplot(bin.loco, aes(x = L3, y = value, colour = L1))
ggloco <- ggloco + geom_line(alpha = 0.5)
ggloco <- ggloco + facet_wrap(~ L2, scales = 'free_y')
ggloco <- ggloco + scale_color_manual(values = cbp)
ggloco <- ggloco + labs(x = 'Time (My)')
ggloco <- ggloco + theme(legend.position = 'none', 
                         axis.title = element_text(size = 19),
                         axis.text = element_text(size = 17))
#ggloco <- ggloco + stat_smooth(method = 'loess', se = FALSE)
ggsave(file = '../doc/figure/loco_bin.png', 
       width = 15, height = 10, plot = ggloco)
