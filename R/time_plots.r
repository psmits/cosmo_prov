library(mgcv)
library(reshape2)
library(ggplot2)
library(scales)

source('../R/cosmo_prov.r')
source('../R/diet_dynamics.r')
source('../R/life_dynamics.r')
source('../R/oxygen_curve.r')

theme_set(theme_bw())
cbp <- c('#E69F00', '#56B4E9', '#009E73', '#F0E442', 
         '#0072B2', '#D55E00', '#CC79A7')

time.match <- cbind(ma = dat$stmid, st = data.frame(dat$stage))
tm <- unique(time.match)

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

# stage
st.dat <- melt(stgraph.bg)
st.dat <- cbind(st.dat, year = tm[match(st.dat$L2, tm[, 2]), 1])
gst <- ggplot(st.dat, aes(x = year, y = value)) + geom_line()
gst <- gst + facet_wrap(~ L1, scales = 'free_y')
gst <- gst + scale_color_manual(values = cbp)
gst <- gst + labs(x = 'Time (My)')
ggsave(file = '../doc/figure/gen_stage.png', 
       width = 15, height = 10, plot = gst)

# standard bin
bin.dat <- melt(win.bg)
bin.dat$L2 <- as.numeric(bin.dat$L2)
ggdat <- ggplot(bin.dat, aes(x = L2, y = value)) + geom_line()
ggdat <- ggdat + facet_wrap(~ L1, scales = 'free_y')
ggdat <- ggdat + scale_color_manual(values = cbp)
#ggdat <- ggdat + theme(legend.position = 'none')
ggdat <- ggdat + labs(x = 'Time (My)')
ggdat <- ggdat + stat_smooth(method = 'loess', se = FALSE)
ggsave(file = '../doc/figure/gen_bin.png', 
       width = 15, height = 10, plot = ggdat)


# diet
# stage
dtst.dat <- melt(stdigr.bg)
dtst.dat <- cbind(dtst.dat, year = tm[match(dtst.dat$L3, tm[, 2]), 1])
dtst.dat$value[is.nan(dtst.dat$value)] <- NA
dtst.dat$value[is.infinite(dtst.dat$value)] <- NA
gdtst <- ggplot(dtst.dat, aes(x = year, y = value, colour = L1))
gdtst <- gdtst + geom_line()
gdtst <- gdtst + facet_wrap(~ L2, scales = 'free_y')
gdtst <- gdtst + scale_color_manual(values = cbp)
gdtst <- gdtst + labs(x = 'Time (My)')
ggsave(file = '../doc/figure/diet_stage.png', 
       width = 15, height = 10, plot = gdtst)

# standard bin
bin.diet <- melt(dtwin.bg)
bin.diet$L3 <- as.numeric(bin.diet$L3)
ggdiet <- ggplot(bin.diet, aes(x = L3, y = value, colour = L1))
ggdiet <- ggdiet + geom_line(alpha = 0.5)
ggdiet <- ggdiet + facet_wrap(~ L2, scales = 'free_y')
ggdiet <- ggdiet + scale_color_manual(values = cbp)
#ggdiet <- ggdiet + theme(legend.position = 'none')
ggdiet <- ggdiet + labs(x = 'Time (My)')
ggdiet <- ggdiet + theme(axis.title = element_text(size = 19),
                         axis.text = element_text(size = 17))
ggdiet <- ggdiet + stat_smooth(method = 'loess', se = FALSE)
ggsave(file = '../doc/figure/diet_bin.png', 
       width = 15, height = 10, plot = ggdiet)


# locomotor
# stage
lfst.dat <- melt(stlfgr.bg)
lfst.dat <- cbind(lfst.dat, year = tm[match(lfst.dat$L3, tm[, 2]), 1])
lfst.dat$value[is.nan(lfst.dat$value)] <- NA
lfst.dat$value[is.infinite(lfst.dat$value)] <- NA
glfst <- ggplot(lfst.dat, aes(x = year, y = value, colour = L1))
glfst <- glfst + geom_line()
glfst <- glfst + facet_wrap(~ L2, scales = 'free_y')
glfst <- glfst + scale_color_manual(values = cbp)
glfst <- glfst + labs(x = 'Time (My)')
ggsave(file = '../doc/figure/loco_stage.png', 
       width = 15, height = 10, plot = glfst)

bin.loco <- melt(lfwin.bg)
bin.loco$L3 <- as.numeric(bin.loco$L3)
ggloco <- ggplot(bin.loco, aes(x = L3, y = value, colour = L1))
ggloco <- ggloco + geom_line(alpha = 0.5)
ggloco <- ggloco + facet_wrap(~ L2, scales = 'free_y')
ggloco <- ggloco + scale_color_manual(values = cbp)
ggloco <- ggloco + labs(x = 'Time (My)')
ggloco <- ggloco + theme(axis.title = element_text(size = 19),
                         axis.text = element_text(size = 17))
#ggloco <- ggloco + stat_smooth(method = 'loess', se = FALSE)
ggsave(file = '../doc/figure/loco_bin.png', 
       width = 15, height = 10, plot = ggloco)
