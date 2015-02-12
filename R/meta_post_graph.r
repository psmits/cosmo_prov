library(ggplot2)
library(reshape2)
library(scales)
library(hexbin)
library(stringr)
library(grid)
library(survival)
library(GGally)
library(plyr)

source('../R/meta_post_sim.r')
nsim <- 100

pairwise.diffs <- function(x) {
  # create column combination pairs
  prs <- cbind(rep(1:ncol(x), each = ncol(x)), 1:ncol(x))
  col.diffs <- prs[prs[, 1] < prs[, 2], , drop = FALSE]
  # do pairwise differences 
  result <- x[, col.diffs[, 1]] - x[, col.diffs[, 2], drop = FALSE]
  # set colnames
  if(is.null(colnames(x)))
    colnames(x) <- 1:ncol(x)
  
  colnames(result) <- paste(colnames(x)[col.diffs[, 1]], ".vs.", 
                            colnames(x)[col.diffs[, 2]], sep = "")
  result
}

theme_set(theme_bw())
cbp <- c('#E69F00', '#56B4E9', '#009E73', '#F0E442', 
         '#0072B2', '#D55E00', '#CC79A7')
theme_update(axis.text = element_text(size = 20),
             axis.title = element_text(size = 30),
             legend.text = element_text(size = 25),
             legend.title = element_text(size = 26),
             legend.key.size = unit(2, 'cm'),
             strip.text = element_text(size = 20))

# discrepency in quantile estimates
# on point
pnt.est <- llply(negbinom.cout, function(y) 
                 aaply(y, .margins = 2, .fun = median))
pnt.est <- cbind(melt(Reduce(rbind, pnt.est)), 
                 label = rep(1:length(pnt.est), 3))[, -1]

mgt.est <- llply(negbinom.cout, function(y) 
                 aaply(y, .margins = 2, .fun = function(x) 
                       quantile(x, c(0.2, 0.8))))
mgt.est <- cbind(Reduce(rbind, mgt.est), 
                 label = rep(1:length(mgt.est), each = 3))
mgt.est <- cbind(data.frame(mgt.est), quant = rownames(mgt.est))
names(mgt.est) <- c('bot', 'top', 'label', 'quant')
mgt.est <- Reduce(rbind, llply(split(mgt.est, mgt.est$quant), 
                               function(x) x[order(x$label), ]))
pnt.est <- cbind(pnt.est, mgt.est[, 1:2])
names(pnt.est) <- c('quant', 'value', 'label', 'bot', 'top')
pnt.est <- pnt.est[!(abs(pnt.est$top - pnt.est$bot) > 100), ]

bad.one <- llply(split(pnt.est, pnt.est$quant), 
                 function(x) which(!(1:33 %in% x$label)))
bad.one <- melt(bad.one)
names(bad.one) <- c('label', 'quant')
bad.one$value <- rep(0, nrow(bad.one))
bad.one$label <- as.numeric(as.character(bad.one$label))

# forward
pnt.for <- llply(negbinom.cfor, function(y)
                 aaply(y, .margins = 2, .fun = median))
pnt.for <- cbind(melt(Reduce(rbind, pnt.for)), 
                 label = rep(1:(length(pnt.for)), 3))[, -1]

mgt.for <- llply(negbinom.cfor, function(y) 
                 aaply(y, .margins = 2, .fun = function(x) 
                       quantile(x, c(0.2, 0.8))))
mgt.for <- cbind(Reduce(rbind, mgt.for), 
                 label = rep(1:(length(mgt.for)), each = 3))
mgt.for <- cbind(data.frame(mgt.for), quant = rownames(mgt.for))
names(mgt.for) <- c('bot', 'top', 'label', 'quant')
mgt.for <- Reduce(rbind, llply(split(mgt.for, mgt.for$quant), 
                               function(x) x[order(x$label), ]))
pnt.for <- cbind(pnt.for, mgt.for[, 1:2])
names(pnt.for) <- c('quant', 'value', 'label', 'bot', 'top')
pnt.for <- pnt.for[!(abs(pnt.for$top - pnt.for$bot) > 100),]

bad.two <- llply(split(pnt.for, pnt.for$quant), 
                 function(x) which(!(1:9 %in% x$label)))
bad.two <- melt(bad.two)
names(bad.two) <- c('label', 'quant')
bad.two$value <- rep(0, nrow(bad.two))
bad.two$label <- as.numeric(as.character(bad.two$label))

pnt.est <- rbind(cbind(pnt.est, mod = rep('est', nrow(pnt.est))),
                 cbind(pnt.for, mod = rep('for', nrow(pnt.for))))

bad.one <- rbind(cbind(bad.one, mod = rep('est', nrow(bad.one))),
                 cbind(bad.two, mod = rep('for', nrow(bad.two))))

disc <- ggplot(pnt.est, aes(x = label, y = value, 
                            ymin = bot, ymax = top, colour = mod))
disc <- disc + geom_hline(aes(yintercept = 0), colour = 'grey', size = 2)
disc <- disc + geom_pointrange(position = position_jitter(height = 0, 
                                                          width = 0.1))
disc <- disc + geom_point(data = bad.one, 
                          mapping = aes(x = label, y = value, 
                                        ymin = value, ymax = value),
                          size = 2, shape = 8, 
                          position = position_jitter(height = 0, 
                                                     width = 0.2))
disc <- disc + facet_grid(quant ~ .)
disc <- disc + labs(x = 'Time (My)', y = 'difference')
disc <- disc + scale_x_reverse()
disc <- disc + scale_colour_manual(values = cbp)



# mass sign/size
eff.mass <- llply(over.coef, function(x) x$beta_mass)
mass.quant <- llply(eff.mass, function(x) quantile(x, c(0.20, 0.5, 0.80)))
mass.quant <- data.frame(Reduce(rbind, mass.quant))
mass.quant <- cbind(mass.quant, label = 2 * (1:(nrow(mass.quant))) - 1)
names(mass.quant) <- c('bot', 'med', 'top', 'label')

efms <- ggplot(mass.quant, aes(x = label, y = med, ymin = bot, ymax = top))
efms <- efms + geom_hline(aes(yintercept = 0), colour = 'grey', size = 2)
efms <- efms + geom_pointrange()
efms <- efms + scale_x_reverse()


# move sign/size
arb.eff <- llply(over.coef, function(x) x$beta_inter)
arb.quant <- llply(arb.eff, function(x) quantile(x, c(0.2, 0.5, 0.8)))
gnd.eff <- llply(over.coef, function(x) x$beta_inter + x$beta_move[, 1])
gnd.quant <- llply(gnd.eff, function(x) quantile(x, c(0.2, 0.5, 0.8)))
scn.eff <- llply(over.coef, function(x) x$beta_inter + x$beta_move[, 2])
scn.quant <- llply(scn.eff, function(x) quantile(x, c(0.2, 0.5, 0.8)))

lamb <- function(x) {
  ll <- length(x)
  temp <- data.frame(Reduce(rbind, x))
  temp <- cbind(temp, label = 2 * (1:ll) + 1)
  names(temp) <- c('bot', 'med', 'top', 'label')
  temp
}
arb.quant <- cbind(lamb(arb.quant), mv = rep('arb', length(arb.quant)))
arb.quant$label <- arb.quant$label + 0.2
gnd.quant <- cbind(lamb(gnd.quant), mv = rep('gnd', length(gnd.quant)))
scn.quant <- cbind(lamb(scn.quant), mv = rep('scn', length(scn.quant)))
scn.quant$label <- scn.quant$label - 0.2
move.eff <- rbind(arb.quant, gnd.quant, scn.quant)

efmv <- ggplot(move.eff, aes(x = label, y = med, 
                             ymin = bot, ymax = top, 
                             colour = mv))
efmv <- efmv + geom_hline(aes(yintercept = 0), colour = 'grey', size = 2)
efmv <- efmv + geom_pointrange()
efmv <- efmv + scale_x_reverse()
efmv <- efmv + scale_colour_manual(values = cbp)


# diet sign/size
crn.eff <- llply(over.coef, function(x) x$beta_inter)
crn.quant <- llply(crn.eff, function(x) quantile(x, c(0.2, 0.5, 0.8)))
hrb.eff <- llply(over.coef, function(x) x$beta_inter + x$beta_diet[, 1])
hrb.quant <- llply(gnd.eff, function(x) quantile(x, c(0.2, 0.5, 0.8)))
ist.eff <- llply(over.coef, function(x) x$beta_inter + x$beta_diet[, 2])
ist.quant <- llply(scn.eff, function(x) quantile(x, c(0.2, 0.5, 0.8)))
omn.eff <- llply(over.coef, function(x) x$beta_inter + x$beta_diet[, 3])
omn.quant <- llply(scn.eff, function(x) quantile(x, c(0.2, 0.5, 0.8)))

crn.quant <- cbind(lamb(crn.quant), mv = rep('crn', length(crn.quant)))
crn.quant$label <- crn.quant$label + 0.4
hrb.quant <- cbind(lamb(hrb.quant), mv = rep('hrb', length(hrb.quant)))
hrb.quant$label <- hrb.quant$label + 0.2
ist.quant <- cbind(lamb(ist.quant), mv = rep('ist', length(ist.quant)))
omn.quant <- cbind(lamb(omn.quant), mv = rep('omn', length(omn.quant)))
omn.quant$label <- omn.quant$label - 0.2

diet.eff <- rbind(crn.quant, hrb.quant, ist.quant, omn.quant)

efdt <- ggplot(diet.eff, aes(x = label, y = med, 
                             ymin = bot, ymax = top, 
                             colour = mv))
efdt <- efdt + geom_hline(aes(yintercept = 0), colour = 'grey', size = 2)
efdt <- efdt + geom_pointrange()
efdt <- efdt + scale_x_reverse()
efdt <- efdt + scale_colour_manual(values = cbp)


# overdispersion in each bin
eff.over <- llply(over.coef, function(x) x$phi)
over.quant <- llply(eff.over, function(x) quantile(x, c(0.20, 0.5, 0.80)))
over.quant <- data.frame(Reduce(rbind, over.quant))
over.quant <- cbind(over.quant, label = 2 * (1:(nrow(over.quant))) - 1)
names(over.quant) <- c('bot', 'med', 'top', 'label')

efov <- ggplot(over.quant, aes(x = label, y = med, ymin = bot, ymax = top))
efov <- efov + geom_hline(aes(yintercept = 1), colour = 'grey', size = 2)
efov <- efov + geom_pointrange()
efov <- efov + scale_x_reverse()
