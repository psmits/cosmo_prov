library(mgcv)
library(ggplot2)
library(reshape2)
library(scales)

source('../R/pa_mung.r')

theme_set(theme_bw())
cbp <- c('#E69F00', '#56B4E9', '#009E73', '#F0E442', 
         '#0072B2', '#D55E00', '#CC79A7')

time.match <- cbind(ma = dat$stmid, st = data.frame(dat$stage))
tm <- unique(time.match)
dat$stage <- factor(dat$stage, levels = tm[order(tm[, 1]), 2])

reldiet <- ddply(dat, .(stage), summarize, 
                 herb = sum(comdiet == 'herb'),
                 omni = sum(comdiet == 'omni'),
                 carni = sum(comdiet == 'carni'))
reldiet <- cbind(ma = tm[, 1], reldiet)
reldiet <- melt(reldiet, id.vars = c('ma', 'stage'))
# plot
rdt <- ggplot(dat, aes(x = stage, fill = comdiet))
rdt <- rdt + geom_bar(position = 'fill')
rdt <- rdt + scale_color_manual(values = cbp)


relloc <- ddply(dat, .(stage), summarize,
                ground = sum(life_habit == 'ground dwelling'),
                fossorial = sum(life_habit == 'fossorial'),
                scansorial = sum(life_habit == 'scansorial'),
                semifoss = sum(life_habit == 'semifossorial'),
                arboreal = sum(life_habit == 'arboreal'),
                staltatorial = sum(life_habit == 'saltatorial'))
# plot
rlf <- ggplot(dat, aes(x = stage, fill = life_habit))
rlf <- rlf + geom_bar(position = 'fill')
rlf <- rlf + scale_color_manual(values = cbp)
