library(mgcv)
library(ggplot2)
library(reshape2)
library(scales)

source('../R/diversity.r')

theme_set(theme_bw())
cbp <- c('#E69F00', '#56B4E9', '#009E73', '#F0E442', 
         '#0072B2', '#D55E00', '#CC79A7')

time.match <- cbind(ma = dat$stmid, st = data.frame(dat$stage))
tm <- unique(time.match)
dat$stage <- factor(dat$stage, levels = tm[order(tm[, 1]), 2])

# relative diet category
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
rdt <- rdt + theme(axis.text.x = element_text(angle = 90))
ggsave(file = '../doc/figure/rel_diet.png', plot = rdt)

# subsampled diet category
subdt <- melt(dietab)
subdt <- cbind(subdt, ma = tm[match(subdt$L2, tm[, 2]), 1])
subdt$L2 <- factor(subdt$L2, levels = tm[order(tm[, 1]), 2])
sdt <- ggplot(subdt, aes(x = ma, y = value, 
                         colour = L1, fill = L1, group = L1))
sdt <- sdt + geom_area(position = 'fill', stat = 'identity')
sdt <- sdt + scale_fill_manual(values = cbp)
sdt <- sdt + scale_color_manual(values = cbp)
ggsave(file = '../doc/figure/sub_diet.png', plot = sdt)

# subsampled at the bin level
subdtbin <- melt(lapply(dtbinab, function(x) {
                        lapply(x, as.numeric)}))
subdtbin$L2 <- as.numeric(subdtbin$L2)
bsdt <- ggplot(subdtbin, aes(x = L2, y = value,
                             colour = L1, fill = L1, group = L1))
bsdt <- bsdt + geom_area(position = 'fill', stat = 'identity')
bsdt <- bsdt + labs(x = 'Time (My)', y = 'relative subsampled richness')
ggsave(file = '../doc/figure/sub_bin_diet.png', plot = sdt)


# relative locomotor category
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
ggsave(file = '../doc/figure/rel_loco.png', plot = rlf)

# subsampled locomotor category
sublf <- melt(locoab)
sublf <- cbind(sublf, ma = tm[match(sublf$L2, tm[, 2]), 1])
sublf$L2 <- factor(sublf$L2, levels = tm[order(tm[, 1]), 2])
slf <- ggplot(sublf, aes(x = ma, y = value, 
                         colour = L1, fill = L1, group = L1))
slf <- slf + geom_area(position = 'fill', stat = 'identity')
slf <- slf + scale_fill_manual(values = cbp)
slf <- slf + scale_color_manual(values = cbp)
ggsave(file = '../doc/figure/sub_loco.png', plot = slf)
