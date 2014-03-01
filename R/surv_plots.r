library(survival)
library(ggplot2)
library(scales)
library(reshape2)
source('../R/step_ribbon.r')

source('../R/surv_analysis.r')

theme_set(theme_bw())
cbp <- c('#E69F00', '#56B4E9', '#009E73', '#F0E442', 
         '#0072B2', '#D55E00', '#CC79A7')

# north america
nakm.curve <- cbind(data.frame(time = na.km$time), surv = na.km$surv)
ggna <- ggplot(nakm.curve, aes(x = time, y = surv))
ggna <- ggna + geom_step()
ggna <- ggna + scale_y_continuous(trans = log_trans())
ggsave(filename = '../doc/figure/km_na.png', plot = ggna,
       height = 10, width = 15)

nn <- c(rep('carni', na.kmd$strata[1]),
        rep('herb', na.kmd$strata[2]),
        rep('insect', na.kmd$strata[3]),
        rep('omni', na.kmd$strata[4]))
nakmd.curve <- cbind(data.frame(time = na.kmd$time), surv = na.kmd$surv, 
                     upper = na.kmd$upper, lower = na.kmd$lower, 
                     diet = nn)
ggnad <- ggplot(nakmd.curve, aes(x = time, y = surv, colour = diet)) 
ggnad <- ggnad + geom_step()
ggnad <- ggnad + geom_ribbon(aes(x = time, ymax = upper, ymin = lower, fill = diet),
                             stat = 'stepribbon', alpha = 0.3, colour = NA)
ggnad <- ggnad + scale_color_manual(values = cbp,
                                    name = 'dietary\ncategory')
ggnad <- ggnad + scale_fill_manual(values = cbp,
                                   name = 'dietary\ncategory')
ggnad <- ggnad + scale_y_continuous(trans = log_trans())
ggnad <- ggnad + coord_cartesian(xlim = c(0, 45))
ggsave(filename = '../doc/figure/km_nad.png', plot = ggnad,
       height = 10, width = 15)

nn <- c(rep('arboreal', na.kml$strata[1]),
        rep('ground dwelling', na.kml$strata[2]),
        rep('scansorial', na.kml$strata[3]))
nakml.curve <- cbind(data.frame(time = na.kml$time), surv = na.kml$surv, 
                     upper = na.kml$upper, lower = na.kml$lower,
                     diet = nn)
ggnal <- ggplot(nakml.curve, aes(x = time, y = surv, colour = diet)) 
ggnal <- ggnal + geom_step()
ggnal <- ggnal + geom_ribbon(aes(x = time, ymax = upper, ymin = lower, fill = diet),
                             stat = 'stepribbon', alpha = 0.3, colour = NA)
ggnal <- ggnal + scale_color_manual(values = cbp,
                                    name = 'locomotor\ncategory')
ggnal <- ggnal + scale_fill_manual(values = cbp,
                                   name = 'locomotor\ncategory')
ggnal <- ggnal + scale_y_continuous(trans = log_trans())
ggnal <- ggnal + coord_cartesian(xlim = c(0, 45))
ggsave(filename = '../doc/figure/km_nal.png', plot = ggnal,
       height = 10, width = 15)


# europe
erkm.curve <- cbind(data.frame(time = er.km$time), surv = er.km$surv)
gger <- ggplot(erkm.curve, aes(x = time, y = surv))
gger <- gger + geom_step()
gger <- gger + scale_y_continuous(trans = log_trans())
ggsave(filename = '../doc/figure/km_er.png', plot = gger,
       height = 10, width = 15)

nn <- c(rep('carni', er.kmd$strata[1]),
        rep('herb', er.kmd$strata[2]),
        rep('insect', er.kmd$strata[3]),
        rep('omni', er.kmd$strata[4]))
erkmd.curve <- cbind(data.frame(time = er.kmd$time), surv = er.kmd$surv, 
                     upper = er.kmd$upper, lower = er.kmd$lower, 
                     diet = nn)
ggerd <- ggplot(erkmd.curve, aes(x = time, y = surv, colour = diet)) 
ggerd <- ggerd + geom_step()
ggerd <- ggerd + geom_ribbon(aes(x = time, ymax = upper, ymin = lower, fill = diet),
                             stat = 'stepribbon', alpha = 0.3, colour = NA)
ggerd <- ggerd + scale_color_manual(values = cbp,
                                    name = 'dietary\ncategory')
ggerd <- ggerd + scale_fill_manual(values = cbp,
                                   name = 'dietary\ncategory')
ggerd <- ggerd + scale_y_continuous(trans = log_trans())
ggerd <- ggerd + coord_cartesian(xlim = c(0, 45))
ggsave(filename = '../doc/figure/km_erd.png', plot = ggerd,
       height = 10, width = 15)

nn <- c(rep('arboreal', er.kml$strata[1]),
        rep('ground dwelling', er.kml$strata[2]),
        rep('scansorial', er.kml$strata[3]))
erkml.curve <- cbind(data.frame(time = er.kml$time), surv = er.kml$surv, 
                     upper = er.kml$upper, lower = er.kml$lower, 
                     diet = nn)
ggerl <- ggplot(erkml.curve, aes(x = time, y = surv, colour = diet)) 
ggerl <- ggerl + geom_step()
ggerl <- ggerl + geom_ribbon(aes(x = time, ymax = upper, ymin = lower, fill = diet),
                             stat = 'stepribbon', alpha = 0.3, colour = NA)
ggerl <- ggerl + scale_color_manual(values = cbp,
                                    name = 'locomotor\ncategory')
ggerl <- ggerl + scale_fill_manual(values = cbp,
                                   name = 'locomotor\ncategory')
ggerl <- ggerl + scale_y_continuous(trans = log_trans())
ggerl <- ggerl + coord_cartesian(xlim = c(0, 45))
ggsave(filename = '../doc/figure/km_erl.png', plot = ggerl,
       height = 10, width = 15)
