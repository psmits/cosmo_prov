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

# generic level
# north america
nag.km.curve <- cbind(data.frame(time = nag.km$time), surv = nag.km$surv)
ggnag <- ggplot(nag.km.curve, aes(x = time, y = surv))
ggnag <- ggnag + geom_step()
ggnag <- ggnag + scale_y_continuous(trans = log_trans())
ggsave(filename = '../doc/figure/km_nag.png', plot = ggnag,
       height = 10, width = 15)
# diet
dd <- c(rep('carni', nag.kmd$strata[1]),
        rep('herb', nag.kmd$strata[2]),
        rep('insect', nag.kmd$strata[3]),
        rep('omni', nag.kmd$strata[4]))
nag.kmd.curve <- cbind(data.frame(time = nag.kmd$time), surv = nag.kmd$surv, 
                     upper = nag.kmd$upper, lower = nag.kmd$lower, 
                     diet = dd)
ggnagd <- ggplot(nag.kmd.curve, aes(x = time, y = surv, colour = diet)) 
ggnagd <- ggnagd + geom_step()
ggnagd <- ggnagd + geom_ribbon(aes(x = time, ymax = upper, ymin = lower, fill = diet),
                             stat = 'stepribbon', alpha = 0.3, colour = NA)
ggnagd <- ggnagd + scale_color_manual(values = cbp,
                                    name = 'dietary\ncategory')
ggnagd <- ggnagd + scale_fill_manual(values = cbp,
                                   name = 'dietary\ncategory')
#ggnagd <- ggnagd + scale_y_continuous(trans = log_trans())
ggnagd <- ggnagd + coord_cartesian(xlim = c(0, 45))
ggsave(filename = '../doc/figure/km_nagd.png', plot = ggnagd,
       height = 10, width = 15)

# locomotor
ll <- c(rep('arboreal', nag.kml$strata[1]),
        rep('ground dwelling', nag.kml$strata[2]),
        rep('scansorial', nag.kml$strata[3]))
nag.kml.curve <- cbind(data.frame(time = nag.kml$time), surv = nag.kml$surv, 
                     upper = nag.kml$upper, lower = nag.kml$lower,
                     loco = ll)
ggnagl <- ggplot(nag.kml.curve, aes(x = time, y = surv, colour = loco)) 
ggnagl <- ggnagl + geom_step()
ggnagl <- ggnagl + geom_ribbon(aes(x = time, ymax = upper, ymin = lower, fill = loco),
                             stat = 'stepribbon', alpha = 0.3, colour = NA)
ggnagl <- ggnagl + scale_color_manual(values = cbp,
                                    name = 'locomotor\ncategory')
ggnagl <- ggnagl + scale_fill_manual(values = cbp,
                                   name = 'locomotor\ncategory')
#ggnagl <- ggnagl + scale_y_continuous(trans = log_trans())
ggnagl <- ggnagl + coord_cartesian(xlim = c(0, 45))
ggsave(filename = '../doc/figure/km_nagl.png', plot = ggnagl,
       height = 10, width = 15)

# diet and locomotor
dl <- c(rep('carnivore, arboreal', nag.kmdl$strata[1]),
        rep('carnivore, ground dwelling', nag.kmdl$strata[2]),
        rep('carnivore, scansorial', nag.kmdl$strata[3]),
        rep('herbivore, arboreal', nag.kmdl$strata[4]),
        rep('herbivore, ground dwelling', nag.kmdl$strata[5]),
        rep('herbivore, scansorial', nag.kmdl$strata[6]),
        rep('insectivore, arboreal', nag.kmdl$strata[7]),
        rep('insectivore, ground dwelling', nag.kmdl$strata[8]),
        rep('insectivore, scansorial', nag.kmdl$strata[9]),
        rep('omnivore, arboreal', nag.kmdl$strata[10]),
        rep('omnivore, ground dwelling', nag.kmdl$strata[11]),
        rep('omnivore, scansorial', nag.kmdl$strata[12]))
nag.kmdl.curve <- cbind(data.frame(time = nag.kmdl$time), surv = nag.kmdl$surv, 
                      upper = nag.kmdl$upper, lower = nag.kmdl$lower,
                      trait = dl)
ggnagdl <- ggplot(nag.kmdl.curve, aes(x = time, y = surv, colour = trait)) 
ggnagdl <- ggnagdl + geom_step()
ggnagdl <- ggnagdl + geom_ribbon(aes(x = time, ymax = upper, ymin = lower, fill = trait),
                              stat = 'stepribbon', alpha = 0.3, colour = NA)
ggnagdl <- ggnagdl + scale_color_manual(values = cbp,
                                     name = 'trait\ncategory')
ggnagdl <- ggnagdl + scale_fill_manual(values = cbp,
                                    name = 'trait\ncategory')
ggnagdl <- ggnagdl + coord_cartesian(xlim = c(0, 45))
ggsave(filename = '../doc/figure/km_nagdl.png', plot = ggnagdl,
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

# generic level
# europe
erg.km.curve <- cbind(data.frame(time = erg.km$time), surv = erg.km$surv)
ggerg <- ggplot(erg.km.curve, aes(x = time, y = surv))
ggerg <- ggerg + geom_step()
ggerg <- ggerg + scale_y_continuous(trans = log_trans())
ggsave(filename = '../doc/figure/km_erg.png', plot = ggerg,
       height = 10, width = 15)
# diet
dd <- c(rep('carni', erg.kmd$strata[1]),
        rep('herb', erg.kmd$strata[2]),
        rep('insect', erg.kmd$strata[3]),
        rep('omni', erg.kmd$strata[4]))
erg.kmd.curve <- cbind(data.frame(time = erg.kmd$time), surv = erg.kmd$surv, 
                     upper = erg.kmd$upper, lower = erg.kmd$lower, 
                     diet = dd)
ggergd <- ggplot(erg.kmd.curve, aes(x = time, y = surv, colour = diet)) 
ggergd <- ggergd + geom_step()
ggergd <- ggergd + geom_ribbon(aes(x = time, ymax = upper, ymin = lower, fill = diet),
                             stat = 'stepribbon', alpha = 0.3, colour = NA)
ggergd <- ggergd + scale_color_manual(values = cbp,
                                    name = 'dietary\ncategory')
ggergd <- ggergd + scale_fill_manual(values = cbp,
                                   name = 'dietary\ncategory')
#ggergd <- ggergd + scale_y_continuous(trans = log_trans())
ggergd <- ggergd + coord_cartesian(xlim = c(0, 45))
ggsave(filename = '../doc/figure/km_ergd.png', plot = ggergd,
       height = 10, width = 15)

# locomotor
ll <- c(rep('arboreal', erg.kml$strata[1]),
        rep('ground dwelling', erg.kml$strata[2]),
        rep('scansorial', erg.kml$strata[3]))
erg.kml.curve <- cbind(data.frame(time = erg.kml$time), surv = erg.kml$surv, 
                     upper = erg.kml$upper, lower = erg.kml$lower,
                     loco = ll)
ggergl <- ggplot(erg.kml.curve, aes(x = time, y = surv, colour = loco)) 
ggergl <- ggergl + geom_step()
ggergl <- ggergl + geom_ribbon(aes(x = time, ymax = upper, ymin = lower, fill = loco),
                             stat = 'stepribbon', alpha = 0.3, colour = NA)
ggergl <- ggergl + scale_color_manual(values = cbp,
                                    name = 'locomotor\ncategory')
ggergl <- ggergl + scale_fill_manual(values = cbp,
                                   name = 'locomotor\ncategory')
#ggergl <- ggergl + scale_y_continuous(trans = log_trans())
ggergl <- ggergl + coord_cartesian(xlim = c(0, 45))
ggsave(filename = '../doc/figure/km_ergl.png', plot = ggergl,
       height = 10, width = 15)

# diet and locomotor
dl <- c(rep('carnivore, arboreal', erg.kmdl$strata[1]),
        rep('carnivore, ground dwelling', erg.kmdl$strata[2]),
        rep('carnivore, scansorial', erg.kmdl$strata[3]),
        rep('herbivore, arboreal', erg.kmdl$strata[4]),
        rep('herbivore, ground dwelling', erg.kmdl$strata[5]),
        rep('herbivore, scansorial', erg.kmdl$strata[6]),
        rep('insectivore, arboreal', erg.kmdl$strata[7]),
        rep('insectivore, ground dwelling', erg.kmdl$strata[8]),
        rep('insectivore, scansorial', erg.kmdl$strata[9]),
        rep('omnivore, arboreal', erg.kmdl$strata[10]),
        rep('omnivore, ground dwelling', erg.kmdl$strata[11]),
        rep('omnivore, scansorial', erg.kmdl$strata[12]))
erg.kmdl.curve <- cbind(data.frame(time = erg.kmdl$time), surv = erg.kmdl$surv, 
                      upper = erg.kmdl$upper, lower = erg.kmdl$lower,
                      trait = dl)
erg.kmdl.curve <- erg.kmdl.curve[-1, ]
erg.kmdl.curve$trait <- factor(as.character(erg.kmdl.curve$trait))
ggergdl <- ggplot(erg.kmdl.curve, aes(x = time, y = surv, colour = trait)) 
ggergdl <- ggergdl + geom_step()
ggergdl <- ggergdl + geom_ribbon(aes(x = time, ymax = upper, ymin = lower, fill = trait),
                              stat = 'stepribbon', alpha = 0.3, colour = NA)
ggergdl <- ggergdl + scale_color_manual(values = cbp,
                                     name = 'trait\ncategory')
ggergdl <- ggergdl + scale_fill_manual(values = cbp,
                                    name = 'trait\ncategory')
ggergdl <- ggergdl + coord_cartesian(xlim = c(0, 45))
ggsave(filename = '../doc/figure/km_ergdl.png', plot = ggergdl,
       height = 10, width = 15)
