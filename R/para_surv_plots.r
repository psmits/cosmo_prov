library(survival)
library(ggplot2)
library(scales)
library(reshape2)
source('../R/step_ribbon.r')

source('../R/surv_parametric.r')

theme_set(theme_bw())
cbp <- c('#A6CEE3', '#B2DF8A', '#FB9a99', '#FDBF6F', '#CAB2D6', '#FFFF99',
         '#1F78B4', '#33A02C', '#E31A1C', '#FF7F00', '#6A3D9A', '#B15928')

# north america
nacurve <- predict(na.exp, 
                   type = 'quantile',
                   p = seq(0.0, 0.99, by = 0.01),
                   se.fit = TRUE)
nacurve <- lapply(nacurve, function(x) x[1, ])
nacurve <- lapply(nacurve, t)
nacurve <- lapply(nacurve, melt)
nacurve <- cbind(fit = nacurve$fit[, -1], se = nacurve$se.fit$value)
nacurve[, 1] <- (100 - nacurve[, 1]) / 100

na <- ggplot(nacurve, aes(x = fit.Var2, y = fit.value))
na <- na + geom_line()
na <- na + geom_ribbon(aes(ymin = fit.value - se, ymax = fit.value + se),
                       alpha = 0.3)
na <- na + coord_flip()
na <- na + labs(y = 'Time', x = 'S(t)')
na <- na + theme(axis.title.y = element_text(angle = 0),
                 axis.text = element_text(size = 20),
                 axis.title = element_text(size = 23),
                 legend.text = element_text(size = 17),
                 legend.title = element_text(size = 19))
ggsave(filename = '../doc/figure/para_na.png', plot = na,
       width = 15, height = 10)

# diet
ndcurve <- predict(nad.exp, newdata = data.frame(diet = c('carni',
                                                          'herb',
                                                          'insect',
                                                          'omni')),
                   type = 'quantile',
                   p = seq(0.0, 0.99, by = 0.01),
                   se.fit = TRUE)
rownames(ndcurve$fit) <- rownames(ndcurve$se.fit) <- c('carni',
                                                       'herb',
                                                       'insect',
                                                       'omni')
ndcurve <- lapply(ndcurve, t)
ndcurve <- lapply(ndcurve, melt)
ndcurve <- cbind(fit = ndcurve$fit, se = ndcurve$se.fit$value)
ndcurve[, 1] <- (100 - ndcurve[, 1]) / 100

nd <- ggplot(ndcurve, aes(x = fit.Var1, y = fit.value, color = fit.Var2))
nd <- nd + geom_line()
nd <- nd + geom_ribbon(aes(ymin = fit.value - se, ymax = fit.value + se,
                       fill = fit.Var2), alpha = 0.3, colour = NA)
nd <- nd + labs(y = 'Time', x = 'S(t)')
nd <- nd + coord_flip()
nd <- nd + scale_fill_manual(values = cbp,
                             name = 'dietary\ncategory')
nd <- nd + scale_colour_manual(values = cbp, 
                               name = 'dietary\ncategory')
nd <- nd + theme(axis.title.y = element_text(angle = 0),
                 axis.text = element_text(size = 20),
                 axis.title = element_text(size = 23),
                 legend.text = element_text(size = 17),
                 legend.title = element_text(size = 19))
ggsave(filename = '../doc/figure/para_na_diet.png', plot = nd,
       width = 15, height = 10)

# locomotor
nlcurve <- predict(nal.exp, newdata = data.frame(move = c('arboreal',
                                                          'ground dwelling',
                                                          'scansorial')),
                   type = 'quantile',
                   p = seq(0.0, 0.99, by = 0.01),
                   se.fit = TRUE)
rownames(nlcurve$fit) <- rownames(nlcurve$se.fit) <- c('arboreal',
                                                       'ground dwelling',
                                                       'scansorial')
nlcurve <- lapply(nlcurve, t)
nlcurve <- lapply(nlcurve, melt)
nlcurve <- cbind(fit = nlcurve$fit, se = nlcurve$se.fit$value)
nlcurve[, 1] <- (100 - nlcurve[, 1]) / 100
nl <- ggplot(nlcurve, aes(x = fit.Var1, y = fit.value, color = fit.Var2))
nl <- nl + geom_line()
nl <- nl + geom_ribbon(aes(ymin = fit.value - se, ymax = fit.value + se,
                       fill = fit.Var2), alpha = 0.3, colour = NA)
nl <- nl + labs(y = 'Time', x = 'S(t)')
nl <- nl + coord_flip()
nl <- nl + scale_fill_manual(values = cbp,
                             name = 'locomotor\ncategory')
nl <- nl + scale_colour_manual(values = cbp, 
                               name = 'locomotor\ncategory')
nl <- nl + theme(axis.title.y = element_text(angle = 0),
                 axis.text = element_text(size = 20),
                 axis.title = element_text(size = 23),
                 legend.text = element_text(size = 17),
                 legend.title = element_text(size = 19))
ggsave(filename = '../doc/figure/para_na_loco.png', plot = nl,
       width = 15, height = 10)

# diet and loco
# this has to be rather iterative
# plot as facet, with each facet as locomotor and the diet spread in each
move <- c('arboreal', 'ground dwelling', 'scansorial')
diet <- c('carni', 'herb', 'insect', 'omni')
out <- list()
for(ii in seq(length(move))) {
  out[[ii]] <- predict(nadl.exp, newdata = data.frame(move = move[ii],
                                                      diet = diet),
                       type = 'quantile',
                       p = seq(0.0, 0.99, by = 0.01),
                       se.fit = TRUE)
}
names(out) <- move
out <- lapply(out, function(x) lapply(x, function(y) {
                                      rownames(y) <- diet
                                      y}))
out <- lapply(out, function(x) lapply(x, t))
out <- lapply(out, function(x) lapply(x, melt))
out <- lapply(out, function(x) cbind(fit = x$fit, se = x$se.fit$value))
out <- cbind(Reduce(rbind, out), 
             move = Reduce(c, Map(function(x, y) rep(x, y), 
                                  names(out), lapply(out, nrow))))
out[, 1] <- (100 - out[, 1]) / 100
names(out) <- c('quant', 'diet', 'surv', 'se', 'move')
ndl <- ggplot(out, aes(x = quant, y = surv, colour = diet))
ndl <- ndl + geom_line()
ndl <- ndl + geom_ribbon(aes(ymin = surv - se, ymax = surv + se, fill = diet),
                         alpha = 0.3, colour = NA)
ndl <- ndl + coord_flip()
ndl <- ndl + labs(y = 'Time', x = 'S(t)')
ndl <- ndl + facet_grid(. ~ move)
ndl <- ndl + theme(axis.title.y = element_text(angle = 0),
                   axis.text = element_text(size = 20),
                   axis.title = element_text(size = 23),
                   legend.text = element_text(size = 17),
                   legend.title = element_text(size = 19))
ggsave(filename = '../doc/figure/para_na_dl.png', plot = ndl,
       width = 15, height = 10)

# europe
ercurve <- predict(er.wei, 
                   type = 'quantile',
                   p = seq(0.0, 0.99, by = 0.01),
                   se.fit = TRUE)
ercurve <- lapply(ercurve, function(x) x[1, ])
ercurve <- lapply(ercurve, t)
ercurve <- lapply(ercurve, melt)
ercurve <- cbind(fit = ercurve$fit[, -1], se = ercurve$se.fit$value)
ercurve[, 1] <- (100 - ercurve[, 1]) / 100

er <- ggplot(ercurve, aes(x = fit.Var2, y = fit.value))
er <- er + geom_line()
er <- er + geom_ribbon(aes(ymin = fit.value - se, ymax = fit.value + se),
                       alpha = 0.3)
er <- er + coord_flip()
er <- er + labs(y = 'Time', x = 'S(t)')
er <- er + theme(axis.title.y = element_text(angle = 0),
                 axis.text = element_text(size = 20),
                 axis.title = element_text(size = 23),
                 legend.text = element_text(size = 17),
                 legend.title = element_text(size = 19))
ggsave(filename = '../doc/figure/para_er.png', plot = er,
       width = 15, height = 10)

# diet
edcurve <- predict(erd.wei, newdata = data.frame(diet = c('carni',
                                                          'herb',
                                                          'insect',
                                                          'omni')),
                   type = 'quantile',
                   p = seq(0.0, 0.99, by = 0.01),
                   se.fit = TRUE)
rownames(edcurve$fit) <- rownames(edcurve$se.fit) <- c('carni',
                                                       'herb',
                                                       'insect',
                                                       'omni')
edcurve <- lapply(edcurve, t)
edcurve <- lapply(edcurve, melt)
edcurve <- cbind(fit = edcurve$fit, se = edcurve$se.fit$value)
edcurve[, 1] <- (100 - edcurve[, 1]) / 100

ed <- ggplot(edcurve, aes(x = fit.Var1, y = fit.value, color = fit.Var2))
ed <- ed + geom_line()
ed <- ed + geom_ribbon(aes(ymin = fit.value - se, ymax = fit.value + se,
                       fill = fit.Var2), alpha = 0.3, colour = NA)
ed <- ed + labs(y = 'Time', x = 'S(t)')
ed <- ed + coord_flip()
ed <- ed + scale_fill_manual(values = cbp,
                             name = 'dietary\ncategory')
ed <- ed + scale_colour_manual(values = cbp, 
                               name = 'dietary\ncategory')
ed <- ed + theme(axis.title.y = element_text(angle = 0),
                 axis.text = element_text(size = 20),
                 axis.title = element_text(size = 23),
                 legend.text = element_text(size = 17),
                 legend.title = element_text(size = 19))
ggsave(filename = '../doc/figure/para_er_diet.png', plot = ed,
       width = 15, height = 10)

# locomotor
elcurve <- predict(erl.wei, newdata = data.frame(move = c('arboreal',
                                                          'ground dwelling',
                                                          'scansorial')),
                   type = 'quantile',
                   p = seq(0.0, 0.99, by = 0.01),
                   se.fit = TRUE)
rownames(elcurve$fit) <- rownames(elcurve$se.fit) <- c('arboreal',
                                                       'ground dwelling',
                                                       'scansorial')
elcurve <- lapply(elcurve, t)
elcurve <- lapply(elcurve, melt)
elcurve <- cbind(fit = elcurve$fit, se = elcurve$se.fit$value)
elcurve[, 1] <- (100 - elcurve[, 1]) / 100
el <- ggplot(elcurve, aes(x = fit.Var1, y = fit.value, color = fit.Var2))
el <- el + geom_line()
el <- el + geom_ribbon(aes(ymin = fit.value - se, ymax = fit.value + se,
                       fill = fit.Var2), alpha = 0.3, colour = NA)
el <- el + labs(y = 'Time', x = 'S(t)')
el <- el + coord_flip()
el <- el + scale_fill_manual(values = cbp,
                             name = 'locomotor\ncategory')
el <- el + scale_colour_manual(values = cbp, 
                               name = 'locomotor\ncategory')
el <- el + theme(axis.title.y = element_text(angle = 0),
                 axis.text = element_text(size = 20),
                 axis.title = element_text(size = 23),
                 legend.text = element_text(size = 17),
                 legend.title = element_text(size = 19))
ggsave(filename = '../doc/figure/para_er_loco.png', plot = el,
       width = 15, height = 10)

# diet and loco
move <- c('arboreal', 'ground dwelling', 'scansorial')
diet <- c('carni', 'herb', 'insect', 'omni')
out <- list()
for(ii in seq(length(move))) {
  out[[ii]] <- predict(erdl.wei, newdata = data.frame(move = move[ii],
                                                      diet = diet),
                       type = 'quantile',
                       p = seq(0.0, 0.99, by = 0.01),
                       se.fit = TRUE)
}
names(out) <- move
out <- lapply(out, function(x) lapply(x, function(y) {
                                      rownames(y) <- diet
                                      y}))
out <- lapply(out, function(x) lapply(x, t))
out <- lapply(out, function(x) lapply(x, melt))
out <- lapply(out, function(x) cbind(fit = x$fit, se = x$se.fit$value))
out <- cbind(Reduce(rbind, out), 
             move = Reduce(c, Map(function(x, y) rep(x, y), 
                                  names(out), lapply(out, nrow))))
out[, 1] <- (100 - out[, 1]) / 100
names(out) <- c('quant', 'diet', 'surv', 'se', 'move')
edl <- ggplot(out, aes(x = quant, y = surv, colour = diet))
edl <- edl + geom_line()
edl <- edl + geom_ribbon(aes(ymin = surv - se, ymax = surv + se, fill = diet),
                         alpha = 0.3, colour = NA)
edl <- edl + coord_flip()
edl <- edl + labs(y = 'Time', x = 'S(t)')
edl <- edl + facet_grid(. ~ move)
edl <- edl + theme(axis.title.y = element_text(angle = 0),
                   axis.text = element_text(size = 20),
                   axis.title = element_text(size = 23),
                   legend.text = element_text(size = 17),
                   legend.title = element_text(size = 19))
ggsave(filename = '../doc/figure/para_er_dl.png', plot = edl,
       width = 15, height = 10)
