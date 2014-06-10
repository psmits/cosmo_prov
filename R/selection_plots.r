library(ggplot2)
library(reshape2)

source('../R/surv_selection.r')

theme_set(theme_bw())
cbp <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
         "#D55E00", "#CC79A7")

# model selection bars
mod.names <- function(mods) {
  out <- lapply(mods, function(x) paste(as.character(x$call$formula[-2], 
                                                     collapse = ' '))[2])
  out <- Map(function(x, y) paste(x, ', ', y$dist,
                                  collapse = ''), 
             out, mods)
  out
}
nwts <- Weights(unlist(lapply(na.mod, AICc)))
nsel <- cbind(data.frame(wts = nwts), nam = unlist(mod.names(na.mod)))
ewts <- Weights(unlist(lapply(er.mod, AICc)))
esel <- cbind(data.frame(wts = ewts), nam = unlist(mod.names(er.mod)))

ngwts <- Weights(unlist(lapply(nagen.mod, AICc)))
ngsel <- cbind(data.frame(wts = ngwts), nam = unlist(mod.names(nagen.mod)))
egwts <- Weights(unlist(lapply(ergen.mod, AICc)))
egsel <- cbind(data.frame(wts = egwts), nam = unlist(mod.names(ergen.mod)))

sel <- rbind(cbind(nsel, loc = rep('NA', nrow(nsel))),
             cbind(esel, loc = rep('Eur', nrow(esel))))
hsel <- rbind(cbind(ngsel, loc = rep('NA', nrow(ngsel))),
              cbind(egsel, loc = rep('Eur', nrow(egsel))))
osel <- rbind(cbind(sel, heir = rep('species', nrow(sel))),
              cbind(hsel, heir = rep('genera', nrow(hsel))))
osel$nam <- factor(osel$nam, levels = unique(as.character(osel$nam)))

gsel <- ggplot(osel, aes(x = loc, y = wts, fill = nam))
gsel <- gsel + geom_bar(stat = 'identity', colour = 'black', width = 0.7)
gsel <- gsel + coord_flip()
gsel <- gsel + facet_grid(heir ~ .)
gsel <- gsel + scale_fill_grey(name = 'Pred., Dist.', start = 0, end = 0.9)
gsel <- gsel + labs(y = 'Akaike Wts', x = '')
gsel <- gsel + theme(axis.title.y = element_text(angle = 0),
                     axis.text = element_text(size = 20),
                     axis.title = element_text(size = 23),
                     legend.text = element_text(size = 10),
                     legend.title = element_text(size = 19),
                     strip.text = element_text(size = 20))
ggsave(filename = '../doc/figure/sel_wts.png', plot = gsel, 
       width = 15, height = 10)


# relative variable importance
# species
imp <- rbind(cbind(na.imp[order(na.imp$pred), ], 
                   base = na.med[order(na.med$pred), 2],
                   loc = rep('NA', nrow(na.imp))), 
             cbind(er.imp[order(er.imp$pred), ], 
                   base = er.med[order(er.med$pred), 2],
                   loc = rep('Eur', nrow(er.imp))))
gimp <- ggplot(imp, aes(x = pred, y = imp, fill = loc))
gimp <- gimp + geom_bar(stat = 'identity', position = 'dodge')
gimp <- gimp + scale_fill_manual(values = cbp, name = 'Region')
gimp <- gimp + geom_errorbar(aes(x = pred, y = base, 
                                 ymax = base, ymin = base),
                             linetype = 'dashed',
                             position = 'dodge',
                             size = 1)
gimp <- gimp + labs(x = 'variable', y = 'rel.\nimp.')
gimp <- gimp + theme(axis.title.y = element_text(angle = 0),
                     axis.text = element_text(size = 20),
                     axis.title = element_text(size = 23),
                     legend.text = element_text(size = 17),
                     legend.title = element_text(size = 19))
ggsave(filename = '../doc/figure/rel_imp.png', plot = gimp, 
       width = 15, height = 10)

# genera
impg <- rbind(cbind(nagen.imp[order(nagen.imp$pred), ], 
                    base = nagen.med[order(nagen.med$pred), 2],
                    loc = rep('NA', nrow(nagen.imp))), 
              cbind(ergen.imp[order(ergen.imp$pred), ], 
                    base = ergen.med[order(ergen.med$pred), 2],
                    loc = rep('Eur', nrow(ergen.imp))))
gimpg <- ggplot(impg, aes(x = pred, y = imp, fill = loc))
gimpg <- gimpg + geom_bar(stat = 'identity', position = 'dodge')
gimpg <- gimpg + scale_fill_manual(values = cbp, name = 'Region')
gimpg <- gimpg + geom_errorbar(aes(x = pred, y = base, 
                                   ymax = base, ymin = base),
                               linetype = 'dashed',
                               position = 'dodge',
                               size = 1)
gimpg <- gimpg + labs(x = 'variable', y = 'rel.\nimp.')
gimpg <- gimpg + theme(axis.title.y = element_text(angle = 0),
                       axis.text = element_text(size = 20),
                       axis.title = element_text(size = 23),
                       legend.text = element_text(size = 17),
                       legend.title = element_text(size = 19))
ggsave(filename = '../doc/figure/rel_imp_gen.png', plot = gimpg, 
       width = 15, height = 10)

# combined species genera
impc <- rbind(cbind(imp, heir = rep('species', nrow(imp))),
              cbind(impg, heir = rep('genera', nrow(impg))))
gimpc <- ggplot(impc, aes(x = pred, y = imp, fill = loc))
gimpc <- gimpc + geom_bar(stat = 'identity', position = 'dodge')
gimpc <- gimpc + scale_fill_manual(values = cbp[-1], name = 'Region')
gimpc <- gimpc + geom_errorbar(aes(x = pred, y = base, 
                                   ymax = base, ymin = base),
                               linetype = 'dashed',
                               position = 'dodge',
                               size = 1)
gimpc <- gimpc + labs(x = 'variable', y = 'rel.\nimp.')
gimpc <- gimpc + facet_grid(. ~ heir, scales = 'free_x')
gimpc <- gimpc + theme(axis.title.y = element_text(angle = 0),
                       axis.text = element_text(size = 20),
                       axis.title = element_text(size = 23),
                       legend.text = element_text(size = 17),
                       legend.title = element_text(size = 19),
                       strip.text = element_text(size = 20))
ggsave(filename = '../doc/figure/rel_imp_com.png', plot = gimpc, 
       width = 15, height = 10)


# model averaging of shape parameter
clean.shape <- function(x) {
  shapes <- laply(x, function(y) y$shape)
  ses <- laply(x, function(y) y$se)
  out <- data.frame(cbind(shapes, ses))
  out
}

na.ks <- clean.shape(na.shape)
er.ks <- clean.shape(er.shape)
ks <- rbind(cbind(na.ks, loc = rep('NA', nrow(na.ks))),
            cbind(er.ks, loc = rep('Eur', nrow(er.ks))))
nag.ks <- clean.shape(nagen.shape)
erg.ks <- clean.shape(ergen.shape)
ksg <- rbind(cbind(nag.ks, loc = rep('NA', nrow(nag.ks))),
             cbind(erg.ks, loc = rep('Eur', nrow(erg.ks))))

ksc <- rbind(cbind(ks, heir = rep('species', nrow(ks))),
             cbind(ksg, heir = rep('genera', nrow(ksg))))

gk <- ggplot(ksc, aes(x = loc, y = shapes))
gk <- gk + geom_pointrange(aes(ymax = shapes + ses, ymin = shapes - ses), 
                           position = position_jitter(w = 0.2, h = 0),
                           size = 1)
gk <- gk + geom_hline(yintercept = 1, size = 1)
gk <- gk + facet_grid(. ~ heir)
gk <- gk + labs(x = '', y = 'shape')
gk <- gk + theme(axis.title.y = element_text(angle = 0),
                 axis.text = element_text(size = 20),
                 axis.title = element_text(size = 23),
                 legend.text = element_text(size = 17),
                 legend.title = element_text(size = 19),
                 strip.text = element_text(size = 20))
ggsave(filename = '../doc/figure/shape_est.png', plot = gk,
       width = 15, height = 10)
