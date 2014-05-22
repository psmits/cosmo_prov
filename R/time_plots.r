library(mgcv)
library(reshape2)
library(ggplot2)
library(scales)

source('../R/na_mung.r')
source('../R/europe_mung.r')

source('../R/cosmo_prov.r')
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

bin.dat$L2[bin.dat$L2 == 'bc'] <- 'BC' 
bin.dat$L2[bin.dat$L2 == 'end'] <- 'E'
bin.dat$L2[bin.dat$L2 == 'avgcoc'] <- 'Occ'
bin.dat$L2[bin.dat$L2 == 'code'] <- 'Code length'

bin.dat$L1[bin.dat$L1 == 'na'] <- 'NA'
bin.dat$L1[bin.dat$L1 == 'eur'] <- 'Europe'

ggdat <- ggplot(bin.dat, aes(x = L3, y = value, linetype = L1)) + geom_line()
ggdat <- ggdat + labs(x = 'Time (My)')
ggdat <- ggdat + scale_linetype_manual(values = c(1,2),
                                       name = 'Region')
ggdat <- ggdat + stat_smooth(method = 'loess', se = FALSE, na.rm = TRUE)
ggdat <- ggdat + facet_wrap(~ L2, scales = 'free_y')
ggdat <- ggdat + theme(axis.title.y = element_text(angle = 0),
                       axis.text = element_text(size = 17),
                       axis.title = element_text(size = 20),
                       legend.text = element_text(size = 15),
                       legend.title = element_text(size = 16),
                       strip.text = element_text(size = 15))
ggsave(file = '../doc/figure/gen_bin.png', 
       width = 15, height = 10, plot = ggdat)


shrink.trait <- function(trait, loc = 'NA', key = 'insect') {
  nat <- lapply(trait, melt)
  ord <- lapply(nat, function(x) {
                apply(x, 2, function(y) any(key %in% y))})
  mat <- Map(function(x, y) {
             names(x)[y] <- 'cat'
             x}, nat, ord)
  mat <- lapply(mat, function(x) {
                oo <- order(x$cat)
                x <- x[oo, ]
                oth <- names(x)[names(x) != 'cat']
                x[, c(oth, 'cat')]})
  mat <- lapply(mat, function(x) {
                names(x) <- c('value', 'age', 'trait')
                x})
  mat <- Map(function(x, y) cbind(x, stat = rep(y, nrow(x))), mat, names(mat))
  mat <- Reduce(rbind, mat)
  mat <- cbind(mat, loc = rep(loc, nrow(mat)))

  mat
}

# diet
dietdf <- list(na = shrink.trait(na.trait$diet),
               eur = shrink.trait(er.trait$diet, 'Eur'))
dietdf <- Reduce(rbind, dietdf)
dietdf$age <- as.numeric(dietdf$age)

ggdiet <- ggplot(dietdf, aes(x = age, y = value, colour = trait))
ggdiet <- ggdiet + geom_line()
ggdiet <- ggdiet + scale_color_manual(values = cbp)
ggdiet <- ggdiet + labs(x = 'Time (My)')
ggdiet <- ggdiet + facet_grid(loc ~ stat, scales = 'free')
ggsave(file = '../doc/figure/diets.png', width = 15, height = 10, plot = ggdiet)


composition <- dietdf[dietdf$stat %in% c('end', 'avgcoc'), ]
composition$stat <- as.character(composition$stat)
composition$stat[composition$stat == 'end'] <- 'E'
composition$stat[composition$stat == 'avgcoc'] <- 'Occ'
composition$loc <- as.character(composition$loc)
composition$loc[composition$loc == 'Eur'] <- 'Europe'
ggcomp <- ggplot(composition, aes(x = age, y = value, colour = trait))
ggcomp <- ggcomp + geom_line()
ggcomp <- ggcomp + scale_color_manual(values = cbp,
                                      name = 'Dietary\nCategory')
ggcomp <- ggcomp + labs(x = 'Time (My)')
ggcomp <- ggcomp + facet_grid(loc ~ stat, scales = 'free')
ggcomp <- ggcomp + theme(axis.title.y = element_text(angle = 0),
                         axis.text = element_text(size = 17),
                         axis.title = element_text(size = 20),
                         legend.text = element_text(size = 15),
                         legend.title = element_text(size = 16),
                         strip.text = element_text(size = 15))
ggsave(file = '../doc/figure/comp_diet.png', width = 15, height = 10, plot = ggcomp)


# just NA
nadt <- shrink.trait(na.trait$diet)
nadt$age <- as.numeric(nadt$age)
nadt$stat <- as.character(nadt$stat)
nadt$stat[nadt$stat == 'bc'] <- 'BC'
nadt$stat[nadt$stat == 'end'] <- 'E'
nadt$stat[nadt$stat == 'avgcoc'] <- 'Occ'
nadt$stat[nadt$stat == 'code'] <- 'Code length'

ggnad <- ggplot(nadt, aes(x = age, y = value, colour = trait))
ggnad <- ggnad + geom_line()
#ggnad <- ggnad + stat_smooth(method = 'loess', se = FALSE, na.rm = TRUE)
ggnad <- ggnad + scale_color_manual(values = cbp,
                                    name = 'Dietary\nCategory')
ggnad <- ggnad + labs(x = 'Time (My)')
ggnad <- ggnad + facet_wrap(~ stat, scales = 'free')
ggnad <- ggnad + theme(axis.title.y = element_text(angle = 0),
                       axis.text = element_text(size = 17),
                       axis.title = element_text(size = 20),
                       legend.text = element_text(size = 15),
                       legend.title = element_text(size = 16),
                       strip.text = element_text(size = 15))
ggsave(file = '../doc/figure/na_dt.png', width = 15, height = 10, plot = ggnad)

# just Eur
erdt <- shrink.trait(er.trait$diet, 'Eur')
erdt$age <- as.numeric(erdt$age)
erdt$stat <- as.character(erdt$stat)
erdt$stat[erdt$stat == 'bc'] <- 'BC'
erdt$stat[erdt$stat == 'end'] <- 'E'
erdt$stat[erdt$stat == 'avgcoc'] <- 'Occ'

ggerd <- ggplot(erdt, aes(x = age, y = value, colour = trait))
ggerd <- ggerd + geom_line()
#ggerd <- ggerd + stat_smooth(method = 'loess', se = FALSE, na.rm = TRUE)
ggerd <- ggerd + scale_color_manual(values = cbp,
                                    name = 'Dietary\nCategory')
ggerd <- ggerd + labs(x = 'Time (My)')
ggerd <- ggerd + facet_wrap(~ stat, scales = 'free')
ggerd <- ggerd + theme(axis.title.y = element_text(angle = 0),
                       axis.text = element_text(size = 17),
                       axis.title = element_text(size = 20),
                       legend.text = element_text(size = 15),
                       legend.title = element_text(size = 16),
                       strip.text = element_text(size = 15))
ggsave(file = '../doc/figure/er_dt.png', width = 15, height = 10, plot = ggerd)


# locomotor
locodf <- list(na = shrink.trait(na.trait$life, key = 'arboreal'),
               eur = shrink.trait(er.trait$life, 'Eur', key = 'arboreal'))
locodf <- Reduce(rbind, locodf)
locodf$age <- as.numeric(locodf$age)

ggloco <- ggplot(locodf, aes(x = age, y = value, colour = trait))
ggloco <- ggloco + geom_line()
ggloco <- ggloco + scale_color_manual(values = cbp)
ggloco <- ggloco + labs(x = 'Time (My)')
ggloco <- ggloco + facet_grid(loc ~ stat, scales = 'free')
ggsave(file = '../doc/figure/locos.png', width = 15, height = 10, plot = ggloco)

grouping <- locodf[locodf$stat %in% c('end', 'avgcoc'), ]
grouping$stat <- as.character(grouping$stat)
grouping$stat[grouping$stat == 'end'] <- 'E'
grouping$stat[grouping$stat == 'avgcoc'] <- 'Occ'
grouping$loc <- as.character(grouping$loc)
grouping$loc[grouping$loc == 'Eur'] <- 'Europe'
ggroup <- ggplot(grouping, aes(x = age, y = value, colour = trait))
ggroup <- ggroup + geom_line()
ggroup <- ggroup + scale_color_manual(values = cbp,
                                      name = 'Locomotor\nCategory')
ggroup <- ggroup + labs(x = 'Time (My)')
ggroup <- ggroup + facet_grid(loc ~ stat, scales = 'free')
ggroup <- ggroup + theme(axis.title.y = element_text(angle = 0),
                         axis.text = element_text(size = 17),
                         axis.title = element_text(size = 20),
                         legend.text = element_text(size = 15),
                         legend.title = element_text(size = 16),
                         strip.text = element_text(size = 15))
ggsave(file = '../doc/figure/comp_loco.png', width = 15, height = 10, plot = ggroup)


# just NA
nalf <- shrink.trait(na.trait$life, key = 'arboreal')
nalf$age <- as.numeric(nalf$age)
nalf$stat <- as.character(nalf$stat)
nalf$stat[nalf$stat == 'bc'] <- 'BC'
nalf$stat[nalf$stat == 'end'] <- 'E'
nalf$stat[nalf$stat == 'avgcoc'] <- 'Occ'
nalf$stat[nalf$stat == 'code'] <- 'Code length'

ggnal <- ggplot(nalf, aes(x = age, y = value, colour = trait))
ggnal <- ggnal + geom_line()
#ggnal <- ggnal + stat_smooth(method = 'loess', se = FALSE, na.rm = TRUE)
ggnal <- ggnal + scale_color_manual(values = cbp,
                                    name = 'Locomotor\nCategory')
ggnal <- ggnal + labs(x = 'Time (My)')
ggnal <- ggnal + facet_wrap(~ stat, scales = 'free')
ggnal <- ggnal + theme(axis.title.y = element_text(angle = 0),
                       axis.text = element_text(size = 17),
                       axis.title = element_text(size = 20),
                       legend.text = element_text(size = 15),
                       legend.title = element_text(size = 16),
                       strip.text = element_text(size = 15))
ggsave(file = '../doc/figure/na_lf.png', width = 15, height = 10, plot = ggnal)

# just Eur
erlf <- shrink.trait(er.trait$life, key = 'arboreal')
erlf$age <- as.numeric(erlf$age)
erlf$stat <- as.character(erlf$stat)
erlf$stat[erlf$stat == 'bc'] <- 'BC'
erlf$stat[erlf$stat == 'end'] <- 'E'
erlf$stat[erlf$stat == 'avgcoc'] <- 'Occ'
erlf$stat[erlf$stat == 'code'] <- 'Code length'

ggerl <- ggplot(erlf, aes(x = age, y = value, colour = trait))
ggerl <- ggerl + geom_line()
#ggerl <- ggerl + stat_smooth(method = 'loess', se = FALSE, na.rm = TRUE)
ggerl <- ggerl + scale_color_manual(values = cbp,
                                    name = 'Locomotor\nCategory')
ggerl <- ggerl + labs(x = 'Time (My)')
ggerl <- ggerl + facet_wrap(~ stat, scales = 'free')
ggerl <- ggerl + theme(axis.title.y = element_text(angle = 0),
                       axis.text = element_text(size = 17),
                       axis.title = element_text(size = 20),
                       legend.text = element_text(size = 15),
                       legend.title = element_text(size = 16),
                       strip.text = element_text(size = 15))
ggsave(file = '../doc/figure/er_lf.png', width = 15, height = 10, plot = ggerl)
