library(survival)
library(ggplot2)
library(scales)
library(reshape2)
source('../R/step_ribbon.r')

source('../R/surv_nonpara.r')

theme_set(theme_bw())
cbp <- c('#A6CEE3', '#B2DF8A', '#FB9a99', '#FDBF6F', '#CAB2D6', '#FFFF99',
         '#1F78B4', '#33A02C', '#E31A1C', '#FF7F00', '#6A3D9A', '#B15928')

# nonparametric K-M survival curves

# regional
# species
nakm <- cbind(data.frame(time = na.species[[1]]$time), 
              surv = na.species[[1]]$surv,
              upper = na.species[[1]]$upper, lower = na.species[[1]]$lower,
              loc = rep('North America', length(na.species[[1]]$time)))
erkm <- cbind(data.frame(time = er.species[[1]]$time), 
              surv = er.species[[1]]$surv,
              upper = er.species[[1]]$upper, lower = er.species[[1]]$lower,
              loc = rep('Europe', length(er.species[[1]]$time)))
species <- rbind(nakm, erkm)
gspec <- ggplot(species, aes(x = time, y = surv, colour = loc))
gspec <- gspec + geom_step(size = 1)
gspec <- gspec + geom_ribbon(aes(x = time, ymax = upper, ymin = lower, fill = loc),
                             stat = 'stepribbon', alpha = 0.3, colour = NA)
gspec <- gspec + scale_color_manual(values = cbp,
                                    name = 'Region')
gspec <- gspec + scale_fill_manual(values = cbp,
                                   name = 'Region')
gspec <- gspec + labs(x = 'Duration (My)', y = 'P(T > t)')
gspec <- gspec + theme(axis.title.y = element_text(angle = 0),
                       axis.text = element_text(size = 20),
                       axis.title = element_text(size = 23),
                       legend.text = element_text(size = 17),
                       legend.title = element_text(size = 19))
#gspec <- gspec + scale_y_continuous(trans = log_trans())
ggsave(filename = '../doc/figure/kms_region.png', plot = gspec,
       height = 10, width = 15)

# genera
nagkm <- cbind(data.frame(time = na.genera[[1]]$time), 
               surv = na.genera[[1]]$surv,
               upper = na.genera[[1]]$upper, lower = na.genera[[1]]$lower,
               loc = rep('North America', length(na.genera[[1]]$time)))
ergkm <- cbind(data.frame(time = er.genera[[1]]$time), 
               surv = er.genera[[1]]$surv,
               upper = er.genera[[1]]$upper, lower = er.genera[[1]]$lower, 
               loc = rep('Europe', length(er.genera[[1]]$time)))
genera <- rbind(nagkm, ergkm)
ggen <- ggplot(genera, aes(x = time, y = surv, colour = loc))
ggen <- ggen + geom_step()
ggen <- ggen + geom_ribbon(aes(x = time, ymax = upper, ymin = lower, fill = loc),
                           stat = 'stepribbon', alpha = 0.3, colour = NA)
ggen <- ggen + scale_color_manual(values = cbp,
                                  name = 'Region')
ggen <- ggen + scale_fill_manual(values = cbp,
                                 name = 'Region')
ggen <- ggen + labs(x = 'Duration (My)', y = 'P(T > t)')
ggen <- ggen + theme(axis.title.y = element_text(angle = 0),
                     axis.text = element_text(size = 20),
                     axis.title = element_text(size = 23),
                     legend.text = element_text(size = 17),
                     legend.title = element_text(size = 19))
#ggen <- ggen + scale_y_continuous(trans = log_trans())
ggsave(filename = '../doc/figure/kmg_region.png', plot = ggen,
       height = 10, width = 15)


# diet
# species
dd <- c(rep('carni', na.species[[2]]$strata[1]),
        rep('herb', na.species[[2]]$strata[2]),
        rep('insect', na.species[[2]]$strata[3]),
        rep('omni', na.species[[2]]$strata[4]))
nakmd <- cbind(data.frame(time = na.species[[2]]$time), 
               surv = na.species[[2]]$surv, 
               upper = na.species[[2]]$upper, 
               lower = na.species[[2]]$lower, 
               diet = dd,
               loc = rep('North America', length(na.species[[2]]$time)))
dd <- c(rep('carni', er.species[[2]]$strata[1]),
        rep('herb', er.species[[2]]$strata[2]),
        rep('insect', er.species[[2]]$strata[3]),
        rep('omni', er.species[[2]]$strata[4]))
erkmd <- cbind(data.frame(time = er.species[[2]]$time), 
               surv = er.species[[2]]$surv, 
               upper = er.species[[2]]$upper, 
               lower = er.species[[2]]$lower, 
               diet = dd,
               loc = rep('Europe', length(er.species[[2]]$time)))
diet <- rbind(nakmd, erkmd)
gdiet <- ggplot(diet, aes(x = time, y = surv, colour = diet))
gdiet <- gdiet + geom_step()
gdiet <- gdiet + geom_ribbon(aes(x = time, ymax = upper, ymin = lower, fill = diet),
                             stat = 'stepribbon', alpha = 0.3, colour = NA)
gdiet <- gdiet + facet_grid(~ loc)
gdiet <- gdiet + scale_color_manual(values = cbp,
                                    name = 'Dietary\nCategory')
gdiet <- gdiet + scale_fill_manual(values = cbp,
                                   name = 'Dietary\nCategory')
gdiet <- gdiet + labs(x = 'Duration (My)', y = 'P(T > t)')
gdiet <- gdiet + theme(axis.title.y = element_text(angle = 0),
                       axis.text = element_text(size = 20),
                       axis.title = element_text(size = 23),
                       legend.text = element_text(size = 17),
                       legend.title = element_text(size = 19))
ggsave(filename = '../doc/figure/diet_sp_region.png', plot = gdiet,
       height = 10, width = 15)

# genera
dd <- c(rep('carni', na.genera[[2]]$strata[1]),
        rep('herb', na.genera[[2]]$strata[2]),
        rep('insect', na.genera[[2]]$strata[3]),
        rep('omni', na.genera[[2]]$strata[4]))
nakmd <- cbind(data.frame(time = na.genera[[2]]$time), 
               surv = na.genera[[2]]$surv, 
               upper = na.genera[[2]]$upper, 
               lower = na.genera[[2]]$lower, 
               diet = dd,
               loc = rep('North America', length(na.genera[[2]]$time)))
dd <- c(rep('carni', er.genera[[2]]$strata[1]),
        rep('herb', er.genera[[2]]$strata[2]),
        rep('insect', er.genera[[2]]$strata[3]),
        rep('omni', er.genera[[2]]$strata[4]))
erkmd <- cbind(data.frame(time = er.genera[[2]]$time), 
               surv = er.genera[[2]]$surv, 
               upper = er.genera[[2]]$upper, 
               lower = er.genera[[2]]$lower, 
               diet = dd,
               loc = rep('Europe', length(er.genera[[2]]$time)))
dietg <- rbind(nakmd, erkmd)
gdietg <- ggplot(dietg, aes(x = time, y = surv, colour = diet))
gdietg <- gdietg + geom_step()
gdietg <- gdietg + geom_ribbon(aes(x = time, ymax = upper, ymin = lower, fill = diet),
                             stat = 'stepribbon', alpha = 0.3, colour = NA)
gdietg <- gdietg + facet_grid(~ loc)
gdietg <- gdietg + scale_color_manual(values = cbp,
                                    name = 'Dietary\nCategory')
gdietg <- gdietg + scale_fill_manual(values = cbp,
                                   name = 'Dietary\nCategory')
gdietg <- gdietg + labs(x = 'Duration (My)', y = 'P(T > t)')
gdietg <- gdietg + theme(axis.title.y = element_text(angle = 0),
                         axis.text = element_text(size = 20),
                         axis.title = element_text(size = 23),
                         legend.text = element_text(size = 17),
                         legend.title = element_text(size = 19))
ggsave(filename = '../doc/figure/diet_gen_region.png', plot = gdietg,
       height = 10, width = 15)


# loco
# species
dd <- c(rep('arboreal', na.species[[3]]$strata[1]),
        rep('ground', na.species[[3]]$strata[2]),
        rep('scansorial', na.species[[3]]$strata[3]))
nakml <- cbind(data.frame(time = na.species[[3]]$time), 
               surv = na.species[[3]]$surv, 
               upper = na.species[[3]]$upper, 
               lower = na.species[[3]]$lower, 
               loco = dd,
               loc = rep('North America', length(na.species[[3]]$time)))
dd <- c(rep('arboreal', er.species[[3]]$strata[1]),
        rep('ground', er.species[[3]]$strata[2]),
        rep('scansorial', er.species[[3]]$strata[3]))
erkml <- cbind(data.frame(time = er.species[[3]]$time), 
               surv = er.species[[3]]$surv, 
               upper = er.species[[3]]$upper, 
               lower = er.species[[3]]$lower, 
               loco = dd,
               loc = rep('Europe', length(er.species[[3]]$time)))
loco <- rbind(nakml, erkml)
gloco <- ggplot(loco, aes(x = time, y = surv, colour = loco))
gloco <- gloco + geom_step()
gloco <- gloco + geom_ribbon(aes(x = time, ymax = upper, ymin = lower, fill = loco),
                             stat = 'stepribbon', alpha = 0.3, colour = NA)
gloco <- gloco + facet_grid(~ loc)
gloco <- gloco + scale_color_manual(values = cbp,
                                    name = 'Locomotor\nCategory')
gloco <- gloco + scale_fill_manual(values = cbp,
                                   name = 'Locomotor\nCategory')
gloco <- gloco + labs(x = 'Duration (My)', y = 'P(T > t)')
gloco <- gloco + theme(axis.title.y = element_text(angle = 0),
                       axis.text = element_text(size = 20),
                       axis.title = element_text(size = 23),
                       legend.text = element_text(size = 17),
                       legend.title = element_text(size = 19))
ggsave(filename = '../doc/figure/loco_sp_region.png', plot = gloco,
       height = 10, width = 15)

# genera
dd <- c(rep('arboreal', na.genera[[3]]$strata[1]),
        rep('ground', na.genera[[3]]$strata[2]),
        rep('scansorial', na.genera[[3]]$strata[3]))
nakml <- cbind(data.frame(time = na.genera[[3]]$time), 
               surv = na.genera[[3]]$surv, 
               upper = na.genera[[3]]$upper, 
               lower = na.genera[[3]]$lower, 
               loco = dd,
               loc = rep('North America', length(na.genera[[3]]$time)))
dd <- c(rep('arboreal', er.genera[[3]]$strata[1]),
        rep('ground', er.genera[[3]]$strata[2]),
        rep('scansorial', er.genera[[3]]$strata[3]))
erkml <- cbind(data.frame(time = er.genera[[3]]$time), 
               surv = er.genera[[3]]$surv, 
               upper = er.genera[[3]]$upper, 
               lower = er.genera[[3]]$lower, 
               loco = dd,
               loc = rep('Europe', length(er.genera[[3]]$time)))
locog <- rbind(nakml, erkml)
glocog <- ggplot(locog, aes(x = time, y = surv, colour = loco))
glocog <- glocog + geom_step()
glocog <- glocog + geom_ribbon(aes(x = time, ymax = upper, ymin = lower, fill = loco),
                             stat = 'stepribbon', alpha = 0.3, colour = NA)
glocog <- glocog + facet_grid(~ loc)
glocog <- glocog + scale_color_manual(values = cbp,
                                    name = 'Locomotor\nCategory')
glocog <- glocog + scale_fill_manual(values = cbp,
                                   name = 'Locomotor\nCategory')
glocog <- glocog + labs(x = 'Duration (My)', y = 'P(T > t)')
glocog <- glocog + theme(axis.title.y = element_text(angle = 0),
                         axis.text = element_text(size = 20),
                         axis.title = element_text(size = 23),
                         legend.text = element_text(size = 17),
                         legend.title = element_text(size = 19))
ggsave(filename = '../doc/figure/loco_gen_region.png', plot = glocog,
       height = 10, width = 15)
