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
nacurve <- predict(na.exp[[1]], 
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
ndcurve <- predict(na.exp[[2]], newdata = data.frame(diet = c('carni',
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
nlcurve <- predict(na.exp[[3]], newdata = data.frame(move = c('arboreal',
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

# generic level
nagcurve <- predict(nag.exp, 
                    type = 'quantile',
                    p = seq(0.0, 0.99, by = 0.01),
                    se.fit = TRUE)
nagcurve <- lapply(nagcurve, function(x) x[1, ])
nagcurve <- lapply(nagcurve, t)
nagcurve <- lapply(nagcurve, melt)
nagcurve <- cbind(fit = nagcurve$fit[, -1], se = nagcurve$se.fit$value)
nagcurve[, 1] <- (100 - nagcurve[, 1]) / 100

nag <- ggplot(nagcurve, aes(x = fit.Var2, y = fit.value))
nag <- nag + geom_line()
nag <- nag + geom_ribbon(aes(ymin = fit.value - se, ymax = fit.value + se),
                         alpha = 0.3)
nag <- nag + coord_flip()
nag <- nag + labs(y = 'Time', x = 'S(t)')
nag <- nag + theme(axis.title.y = element_text(angle = 0),
                   axis.text = element_text(size = 20),
                   axis.title = element_text(size = 23),
                   legend.text = element_text(size = 17),
                   legend.title = element_text(size = 19))
ggsave(filename = '../doc/figure/para_nag.png', plot = nag,
       width = 15, height = 10)

# diet
ngdcurve <- predict(nagd.exp, newdata = data.frame(diet = c('carni',
                                                            'herb',
                                                            'insect',
                                                            'omni')),
                    type = 'quantile',
                    p = seq(0.0, 0.99, by = 0.01),
                    se.fit = TRUE)
rownames(ngdcurve$fit) <- rownames(ngdcurve$se.fit) <- c('carni',
                                                         'herb',
                                                         'insect',
                                                         'omni')
ngdcurve <- lapply(ngdcurve, t)
ngdcurve <- lapply(ngdcurve, melt)
ngdcurve <- cbind(fit = ngdcurve$fit, se = ngdcurve$se.fit$value)
ngdcurve[, 1] <- (100 - ngdcurve[, 1]) / 100

ngd <- ggplot(ngdcurve, aes(x = fit.Var1, y = fit.value, color = fit.Var2))
ngd <- ngd + geom_line()
ngd <- ngd + geom_ribbon(aes(ymin = fit.value - se, ymax = fit.value + se,
                             fill = fit.Var2), alpha = 0.3, colour = NA)
ngd <- ngd + labs(y = 'Time', x = 'S(t)')
ngd <- ngd + coord_flip()
ngd <- ngd + scale_fill_manual(values = cbp,
                               name = 'dietary\ncategory')
ngd <- ngd + scale_colour_manual(values = cbp, 
                                 name = 'dietary\ncategory')
ngd <- ngd + theme(axis.title.y = element_text(angle = 0),
                   axis.text = element_text(size = 20),
                   axis.title = element_text(size = 23),
                   legend.text = element_text(size = 17),
                   legend.title = element_text(size = 19))
ggsave(filename = '../doc/figure/para_nag_diet.png', plot = ngd,
       width = 15, height = 10)

# locomotor
nglcurve <- predict(nagl.exp, newdata = data.frame(move = c('arboreal',
                                                            'ground dwelling',
                                                            'scansorial')),
                    type = 'quantile',
                    p = seq(0.0, 0.99, by = 0.01),
                    se.fit = TRUE)
rownames(nglcurve$fit) <- rownames(nglcurve$se.fit) <- c('arboreal',
                                                         'ground dwelling',
                                                         'scansorial')
nglcurve <- lapply(nglcurve, t)
nglcurve <- lapply(nglcurve, melt)
nglcurve <- cbind(fit = nglcurve$fit, se = nglcurve$se.fit$value)
nglcurve[, 1] <- (100 - nglcurve[, 1]) / 100
ngl <- ggplot(nglcurve, aes(x = fit.Var1, y = fit.value, color = fit.Var2))
ngl <- ngl + geom_line()
ngl <- ngl + geom_ribbon(aes(ymin = fit.value - se, ymax = fit.value + se,
                             fill = fit.Var2), alpha = 0.3, colour = NA)
ngl <- ngl + labs(y = 'Time', x = 'S(t)')
ngl <- ngl + coord_flip()
ngl <- ngl + scale_fill_manual(values = cbp,
                               name = 'locomotor\ncategory')
ngl <- ngl + scale_colour_manual(values = cbp, 
                                 name = 'locomotor\ncategory')
ngl <- ngl + theme(axis.title.y = element_text(angle = 0),
                   axis.text = element_text(size = 20),
                   axis.title = element_text(size = 23),
                   legend.text = element_text(size = 17),
                   legend.title = element_text(size = 19))
ggsave(filename = '../doc/figure/para_nal_loco.png', plot = ngl,
       width = 15, height = 10)

# diet and loco
# this has to be rather iterative
# plot as facet, with each facet as locomotor and the diet spread in each
move <- c('arboreal', 'ground dwelling', 'scansorial')
diet <- c('carni', 'herb', 'insect', 'omni')
out <- list()
for(ii in seq(length(move))) {
  out[[ii]] <- predict(nagdl.exp, newdata = data.frame(move = move[ii],
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
ngdl <- ggplot(out, aes(x = quant, y = surv, colour = diet))
ngdl <- ngdl + geom_line()
ngdl <- ngdl + geom_ribbon(aes(ymin = surv - se, ymax = surv + se, fill = diet),
                           alpha = 0.3, colour = NA)
ngdl <- ngdl + coord_flip()
ngdl <- ngdl + labs(y = 'Time', x = 'S(t)')
ngdl <- ngdl + facet_grid(. ~ move)
ngdl <- ngdl + theme(axis.title.y = element_text(angle = 0),
                     axis.text = element_text(size = 20),
                     axis.title = element_text(size = 23),
                     legend.text = element_text(size = 17),
                     legend.title = element_text(size = 19))
ggsave(filename = '../doc/figure/para_nag_dl.png', plot = ngdl,
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

# generic level
ergcurve <- predict(erg.exp, 
                    type = 'quantile',
                    p = seq(0.0, 0.99, by = 0.01),
                    se.fit = TRUE)
ergcurve <- lapply(ergcurve, function(x) x[1, ])
ergcurve <- lapply(ergcurve, t)
ergcurve <- lapply(ergcurve, melt)
ergcurve <- cbind(fit = ergcurve$fit[, -1], se = ergcurve$se.fit$value)
ergcurve[, 1] <- (100 - ergcurve[, 1]) / 100

erg <- ggplot(ergcurve, aes(x = fit.Var2, y = fit.value))
erg <- erg + geom_line()
erg <- erg + geom_ribbon(aes(ymin = fit.value - se, ymax = fit.value + se),
                         alpha = 0.3)
erg <- erg + coord_flip()
erg <- erg + labs(y = 'Time', x = 'S(t)')
erg <- erg + theme(axis.title.y = element_text(angle = 0),
                   axis.text = element_text(size = 20),
                   axis.title = element_text(size = 23),
                   legend.text = element_text(size = 17),
                   legend.title = element_text(size = 19))
ggsave(filename = '../doc/figure/para_erg.png', plot = erg,
       width = 15, height = 10)

# diet
egdcurve <- predict(ergd.exp, newdata = data.frame(diet = c('carni',
                                                            'herb',
                                                            'insect',
                                                            'omni')),
                    type = 'quantile',
                    p = seq(0.0, 0.99, by = 0.01),
                    se.fit = TRUE)
rownames(egdcurve$fit) <- rownames(egdcurve$se.fit) <- c('carni',
                                                         'herb',
                                                         'insect',
                                                         'omni')
egdcurve <- lapply(egdcurve, t)
egdcurve <- lapply(egdcurve, melt)
egdcurve <- cbind(fit = egdcurve$fit, se = egdcurve$se.fit$value)
egdcurve[, 1] <- (100 - egdcurve[, 1]) / 100

egd <- ggplot(egdcurve, aes(x = fit.Var1, y = fit.value, color = fit.Var2))
egd <- egd + geom_line()
egd <- egd + geom_ribbon(aes(ymin = fit.value - se, ymax = fit.value + se,
                             fill = fit.Var2), alpha = 0.3, colour = NA)
egd <- egd + labs(y = 'Time', x = 'S(t)')
egd <- egd + coord_flip()
egd <- egd + scale_fill_manual(values = cbp,
                               name = 'dietary\ncategory')
egd <- egd + scale_colour_manual(values = cbp, 
                                 name = 'dietary\ncategory')
egd <- egd + theme(axis.title.y = element_text(angle = 0),
                   axis.text = element_text(size = 20),
                   axis.title = element_text(size = 23),
                   legend.text = element_text(size = 17),
                   legend.title = element_text(size = 19))
ggsave(filename = '../doc/figure/para_erg_diet.png', plot = egd,
       width = 15, height = 10)

# locomotor
eglcurve <- predict(ergl.exp, newdata = data.frame(move = c('arboreal',
                                                            'ground dwelling',
                                                            'scansorial')),
                    type = 'quantile',
                    p = seq(0.0, 0.99, by = 0.01),
                    se.fit = TRUE)
rownames(eglcurve$fit) <- rownames(eglcurve$se.fit) <- c('arboreal',
                                                         'ground dwelling',
                                                         'scansorial')
eglcurve <- lapply(eglcurve, t)
eglcurve <- lapply(eglcurve, melt)
eglcurve <- cbind(fit = eglcurve$fit, se = eglcurve$se.fit$value)
eglcurve[, 1] <- (100 - eglcurve[, 1]) / 100
egl <- ggplot(eglcurve, aes(x = fit.Var1, y = fit.value, color = fit.Var2))
egl <- egl + geom_line()
egl <- egl + geom_ribbon(aes(ymin = fit.value - se, ymax = fit.value + se,
                             fill = fit.Var2), alpha = 0.3, colour = NA)
egl <- egl + labs(y = 'Time', x = 'S(t)')
egl <- egl + coord_flip()
egl <- egl + scale_fill_manual(values = cbp,
                               name = 'locomotor\ncategory')
egl <- egl + scale_colour_manual(values = cbp, 
                                 name = 'locomotor\ncategory')
egl <- egl + theme(axis.title.y = element_text(angle = 0),
                   axis.text = element_text(size = 20),
                   axis.title = element_text(size = 23),
                   legend.text = element_text(size = 17),
                   legend.title = element_text(size = 19))
ggsave(filename = '../doc/figure/para_nal_loco.png', plot = egl,
       width = 15, height = 10)

# diet and loco
# this has to be rather iterative
# plot as facet, with each facet as locomotor and the diet spread in each
move <- c('arboreal', 'ground dwelling', 'scansorial')
diet <- c('carni', 'herb', 'insect', 'omni')
out <- list()
for(ii in seq(length(move))) {
  out[[ii]] <- predict(ergdl.exp, newdata = data.frame(move = move[ii],
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
egdl <- ggplot(out, aes(x = quant, y = surv, colour = diet))
egdl <- egdl + geom_line()
egdl <- egdl + geom_ribbon(aes(ymin = surv - se, ymax = surv + se, fill = diet),
                           alpha = 0.3, colour = NA)
egdl <- egdl + coord_flip()
egdl <- egdl + labs(y = 'Time', x = 'S(t)')
egdl <- egdl + facet_grid(. ~ move)
egdl <- egdl + theme(axis.title.y = element_text(angle = 0),
                     axis.text = element_text(size = 20),
                     axis.title = element_text(size = 23),
                     legend.text = element_text(size = 17),
                     legend.title = element_text(size = 19))
ggsave(filename = '../doc/figure/para_erg_dl.png', plot = egdl,
       width = 15, height = 10)
