library(survival)
library(ggplot2)
library(scales)
library(reshape2)
source('../R/step_ribbon.r')

source('../R/surv_nonpara.r')

theme_set(theme_bw())
cbp <- c('#A6CEE3', '#B2DF8A', '#FB9a99', '#FDBF6F', '#CAB2D6', '#FFFF99',
         '#1F78B4', '#33A02C', '#E31A1C', '#FF7F00', '#6A3D9A', '#B15928')

# north america
nakm.curve <- cbind(data.frame(time = na.km$time), surv = na.km$surv)
ggna <- ggplot(nakm.curve, aes(x = time, y = surv))
ggna <- ggna + geom_step()
ggna <- ggna + scale_y_continuous(trans = log_trans())
ggsave(filename = '../doc/figure/km_na.png', plot = ggna,
       height = 10, width = 15)
# diet
dd <- c(rep('carni', na.kmd$strata[1]),
        rep('herb', na.kmd$strata[2]),
        rep('insect', na.kmd$strata[3]),
        rep('omni', na.kmd$strata[4]))
nakmd.curve <- cbind(data.frame(time = na.kmd$time), surv = na.kmd$surv, 
                     upper = na.kmd$upper, lower = na.kmd$lower, 
                     diet = dd)
ggnad <- ggplot(nakmd.curve, aes(x = time, y = surv, colour = diet)) 
ggnad <- ggnad + geom_step()
ggnad <- ggnad + geom_ribbon(aes(x = time, ymax = upper, ymin = lower, fill = diet),
                             stat = 'stepribbon', alpha = 0.3, colour = NA)
ggnad <- ggnad + scale_color_manual(values = cbp,
                                    name = 'dietary\ncategory')
ggnad <- ggnad + scale_fill_manual(values = cbp,
                                   name = 'dietary\ncategory')
#ggnad <- ggnad + scale_y_continuous(trans = log_trans())
ggnad <- ggnad + coord_cartesian(xlim = c(0, 45))
ggsave(filename = '../doc/figure/km_nad.png', plot = ggnad,
       height = 10, width = 15)

# locomotor
ll <- c(rep('arboreal', na.kml$strata[1]),
        rep('ground dwelling', na.kml$strata[2]),
        rep('scansorial', na.kml$strata[3]))
nakml.curve <- cbind(data.frame(time = na.kml$time), surv = na.kml$surv, 
                     upper = na.kml$upper, lower = na.kml$lower,
                     loco = ll)
ggnal <- ggplot(nakml.curve, aes(x = time, y = surv, colour = loco)) 
ggnal <- ggnal + geom_step()
ggnal <- ggnal + geom_ribbon(aes(x = time, ymax = upper, ymin = lower, fill = loco),
                             stat = 'stepribbon', alpha = 0.3, colour = NA)
ggnal <- ggnal + scale_color_manual(values = cbp,
                                    name = 'locomotor\ncategory')
ggnal <- ggnal + scale_fill_manual(values = cbp,
                                   name = 'locomotor\ncategory')
#ggnal <- ggnal + scale_y_continuous(trans = log_trans())
ggnal <- ggnal + coord_cartesian(xlim = c(0, 45))
ggsave(filename = '../doc/figure/km_nal.png', plot = ggnal,
       height = 10, width = 15)

# diet and locomotor
dl <- c(rep('carnivore, arboreal', na.kmdl$strata[1]),
        rep('carnivore, ground dwelling', na.kmdl$strata[2]),
        rep('carnivore, scansorial', na.kmdl$strata[3]),
        rep('herbivore, arboreal', na.kmdl$strata[4]),
        rep('herbivore, ground dwelling', na.kmdl$strata[5]),
        rep('herbivore, scansorial', na.kmdl$strata[6]),
        rep('insectivore, arboreal', na.kmdl$strata[7]),
        rep('insectivore, ground dwelling', na.kmdl$strata[8]),
        rep('insectivore, scansorial', na.kmdl$strata[9]),
        rep('omnivore, arboreal', na.kmdl$strata[10]),
        rep('omnivore, ground dwelling', na.kmdl$strata[11]),
        rep('omnivore, scansorial', na.kmdl$strata[12]))
nakmdl.curve <- cbind(data.frame(time = na.kmdl$time), surv = na.kmdl$surv, 
                      upper = na.kmdl$upper, lower = na.kmdl$lower,
                      trait = dl)
ggnadl <- ggplot(nakmdl.curve, aes(x = time, y = surv, colour = trait)) 
ggnadl <- ggnadl + geom_step()
ggnadl <- ggnadl + geom_ribbon(aes(x = time, ymax = upper, ymin = lower, fill = trait),
                              stat = 'stepribbon', alpha = 0.3, colour = NA)
ggnadl <- ggnadl + scale_color_manual(values = cbp,
                                     name = 'trait\ncategory')
ggnadl <- ggnadl + scale_fill_manual(values = cbp,
                                    name = 'trait\ncategory')
ggnadl <- ggnadl + coord_cartesian(xlim = c(0, 45))
ggsave(filename = '../doc/figure/km_nadl.png', plot = ggnadl,
       height = 10, width = 15)


# europe
erkm.curve <- cbind(data.frame(time = er.km$time), surv = er.km$surv)
gger <- ggplot(erkm.curve, aes(x = time, y = surv))
gger <- gger + geom_step()
gger <- gger + scale_y_continuous(trans = log_trans())
ggsave(filename = '../doc/figure/km_er.png', plot = gger,
       height = 10, width = 15)

# diet
dd <- c(rep('carni', er.kmd$strata[1]),
        rep('herb', er.kmd$strata[2]),
        rep('insect', er.kmd$strata[3]),
        rep('omni', er.kmd$strata[4]))
erkmd.curve <- cbind(data.frame(time = er.kmd$time), surv = er.kmd$surv, 
                     upper = er.kmd$upper, lower = er.kmd$lower, 
                     diet = dd)
ggerd <- ggplot(erkmd.curve, aes(x = time, y = surv, colour = diet)) 
ggerd <- ggerd + geom_step()
ggerd <- ggerd + geom_ribbon(aes(x = time, ymax = upper, ymin = lower, fill = diet),
                             stat = 'stepribbon', alpha = 0.3, colour = NA)
ggerd <- ggerd + scale_color_manual(values = cbp,
                                    name = 'dietary\ncategory')
ggerd <- ggerd + scale_fill_manual(values = cbp,
                                   name = 'dietary\ncategory')
ggerd <- ggerd + coord_cartesian(xlim = c(0, 45))
ggsave(filename = '../doc/figure/km_erd.png', plot = ggerd,
       height = 10, width = 15)

# locomotor
ll <- c(rep('arboreal', er.kml$strata[1]),
        rep('ground dwelling', er.kml$strata[2]),
        rep('scansorial', er.kml$strata[3]))
erkml.curve <- cbind(data.frame(time = er.kml$time), surv = er.kml$surv, 
                     upper = er.kml$upper, lower = er.kml$lower, 
                     loco = ll)
ggerl <- ggplot(erkml.curve, aes(x = time, y = surv, colour = loco)) 
ggerl <- ggerl + geom_step()
ggerl <- ggerl + geom_ribbon(aes(x = time, ymax = upper, ymin = lower, fill = loco),
                             stat = 'stepribbon', alpha = 0.3, colour = NA)
ggerl <- ggerl + scale_color_manual(values = cbp,
                                    name = 'locomotor\ncategory')
ggerl <- ggerl + scale_fill_manual(values = cbp,
                                   name = 'locomotor\ncategory')
ggerl <- ggerl + coord_cartesian(xlim = c(0, 45))
ggsave(filename = '../doc/figure/km_erl.png', plot = ggerl,
       height = 10, width = 15)

# diet and locomotor
dl <- c(rep('carnivore, arboreal', er.kmdl$strata[1]),
        rep('carnivore, ground dwelling', er.kmdl$strata[2]),
        rep('carnivore, scansorial', er.kmdl$strata[3]),
        rep('herbivore, arboreal', er.kmdl$strata[4]),
        rep('herbivore, ground dwelling', er.kmdl$strata[5]),
        rep('herbivore, scansorial', er.kmdl$strata[6]),
        rep('insectivore, arboreal', er.kmdl$strata[7]),
        rep('insectivore, ground dwelling', er.kmdl$strata[8]),
        rep('insectivore, scansorial', er.kmdl$strata[9]),
        rep('omnivore, arboreal', er.kmdl$strata[10]),
        rep('omnivore, ground dwelling', er.kmdl$strata[11]),
        rep('omnivore, scansorial', er.kmdl$strata[12]))
erkmdl.curve <- cbind(data.frame(time = er.kmdl$time), surv = er.kmdl$surv, 
                      upper = er.kmdl$upper, lower = er.kmdl$lower,
                      trait = dl)
ggerdl <- ggplot(erkmdl.curve, aes(x = time, y = surv, colour = trait)) 
ggerdl <- ggerdl + geom_step()
ggerdl <- ggerdl + geom_ribbon(aes(x = time, ymax = upper, ymin = lower, fill = trait),
                              stat = 'stepribbon', alpha = 0.3, colour = NA)
ggerdl <- ggerdl + scale_color_manual(values = cbp,
                                     name = 'trait\ncategory')
ggerdl <- ggerdl + scale_fill_manual(values = cbp,
                                    name = 'trait\ncategory')
ggerdl <- ggerdl + coord_cartesian(xlim = c(0, 45))
ggsave(filename = '../doc/figure/km_erdl.png', plot = ggerdl,
       height = 10, width = 15)
