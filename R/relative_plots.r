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
reldiet <- cbind(ma = tm[match(reldiet$stage, tm[, 2]), 1], reldiet)
reldiet <- melt(reldiet, id.vars = c('ma', 'stage'))
# plot
rdt <- ggplot(reldiet, aes(x = ma, y = value,
                             colour = variable, 
                             fill = variable, 
                             group = variable))
rdt <- rdt + geom_area(position = 'fill', stat = 'identity')
rdt <- rdt + scale_color_manual(values = cbp)
rdt <- rdt + scale_fill_manual(values = cbp)
ggsave(file = '../doc/figure/rel_diet.png', plot = rdt)

# subsampled diet category
subdt <- melt(dietab)
subdt <- cbind(subdt, ma = tm[match(subdt$L2, tm[, 2]), 1])
subdt$L2 <- factor(subdt$L2, levels = tm[order(tm[, 1]), 2])
subdt$L1 <- factor(subdt$L1, levels = levels(reldiet$variable))
sdt <- ggplot(subdt, aes(x = ma, y = value, 
                         colour = L1, fill = L1, group = L1))
sdt <- sdt + geom_area(position = 'fill', stat = 'identity')
sdt <- sdt + scale_fill_manual(values = cbp)
sdt <- sdt + scale_color_manual(values = cbp)
ggsave(file = '../doc/figure/sub_diet.png', plot = sdt)

# faceted version of the the stage plot
rrsub <- subdt[, c(4, 2, 3, 1)]
names(rrsub) <- names(reldiet)
relcombo <- rbind(cbind(reldiet, cl = rep('raw', nrow(reldiet))),
                  cbind(rrsub, cl = rep('sub', nrow(subdt))))
relcombo <- relcombo[order(relcombo$variable), ]
relcb <- ggplot(relcombo, aes(x = ma, y = value,
                              colour = variable, 
                              fill = variable, 
                              group = variable))
relcb <- relcb + geom_area(position = 'fill', stat = 'identity')
relcb <- relcb + scale_fill_manual(values = cbp)
relcb <- relcb + scale_color_manual(values = cbp)
relcb <- relcb + facet_wrap(~ cl, scales = 'free')
ggsave(file = '../doc/figure/facet_diet.png', plot = relcb)


# bin level
relbin <- ddply(dat, .(bins), summarize,
                herb = sum(comdiet == 'herb'),
                omni = sum(comdiet == 'omni'),
                carni = sum(comdiet == 'carni'))
mrb <- melt(relbin, id.vars = 'bins')
mrb <- mrb[, c(3, 1, 2)]

# subsampled at the bin level
subdtbin <- melt(lapply(dtbinab, function(x) {
                        lapply(x, as.numeric)}))
subdtbin$L2 <- as.numeric(subdtbin$L2)
bsdt <- ggplot(subdtbin, aes(x = L2, y = value,
                             colour = L1, fill = L1, group = L1))
bsdt <- bsdt + geom_area(position = 'fill', stat = 'identity')
bsdt <- bsdt + labs(x = 'Time (My)', y = 'relative subsampled richness')
ggsave(file = '../doc/figure/sub_bin_diet.png', plot = sdt)

names(subdtbin) <- names(mrb)
relmix <- rbind(cbind(mrb, cl = rep('raw', nrow(mrb))),
                cbind(subdtbin, cl = rep('sub', nrow(subdtbin))))
relmix$variable <- factor(relmix$variable, 
                          levels = c('carni', 'herb', 'omni'))
relmix <- relmix[order(relmix$variable), ]
#relmix <- relmix[order(relmix$variable), ]
#relmix <- relmix[relmix$variable != 'omni', ]
relmix <- Reduce(rbind, 
                 lapply(split(relmix, relmix$cl), function(x) {
                        rms <- x[is.na(x[, 1]), 2]
                        x[!(x[, 2] %in% c(rms, 6, 8)), ]
                             })
                 )
mix <- ggplot(relmix, aes(x = bins, y = value,
                          colour = variable, 
                          fill = variable,
                          group = variable))
mix <- mix + geom_area(position = 'fill', stat = 'identity')
mix <- mix + labs(x = 'Time (My)', y = 'relative richness')
mix <- mix + scale_fill_manual(values = cbp)
mix <- mix + scale_color_manual(values = cbp)
mix <- mix + facet_wrap(~ cl, scales = 'free')
ggsave(file = '../doc/figure/facet_mix.png', 
       width = 15, height = 10, plot = mix)



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
