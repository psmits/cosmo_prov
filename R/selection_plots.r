library(ggplot2)
library(reshape2)

source('../R/surv_selection.r')

theme_set(theme_bw())
cbp <- c('#A6CEE3', '#B2DF8A', '#FB9a99', '#FDBF6F', '#CAB2D6', '#FFFF99',
         '#1F78B4', '#33A02C', '#E31A1C', '#FF7F00', '#6A3D9A', '#B15928')

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
gimpc <- gimpc + scale_fill_manual(values = cbp, name = 'Region')
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
