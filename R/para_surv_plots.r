library(survival)
library(ggplot2)
library(scales)
library(reshape2)
source('../R/step_ribbon.r')

source('../R/surv_parametric.r')

theme_set(theme_bw())
cbp <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
         "#D55E00", "#CC79A7")

# north america and europe
# species
nacurve <- predict(na.exp[[1]], 
                   type = 'quantile',
                   p = seq(0.01, 0.99, by = 0.01),
                   se.fit = TRUE)
nacurve <- lapply(nacurve, function(x) x[1, ])
nacurve <- lapply(nacurve, t)
nacurve <- lapply(nacurve, melt)
nacurve <- cbind(fit = nacurve$fit[, -1], se = nacurve$se.fit$value)
nacurve[, 1] <- (100 - nacurve[, 1]) / 100

ercurve <- predict(er.wei[[1]], 
                   type = 'quantile',
                   p = seq(0.01, 0.99, by = 0.01),
                   se.fit = TRUE)
ercurve <- lapply(ercurve, function(x) x[1, ])
ercurve <- lapply(ercurve, t)
ercurve <- lapply(ercurve, melt)
ercurve <- cbind(fit = ercurve$fit[, -1], se = ercurve$se.fit$value)
ercurve[, 1] <- (100 - ercurve[, 1]) / 100

# genus
nagcurve <- predict(nagen.exp[[1]], 
                    type = 'quantile',
                    p = seq(0.01, 0.99, by = 0.01),
                    se.fit = TRUE)
nagcurve <- lapply(nagcurve, function(x) x[1, ])
nagcurve <- lapply(nagcurve, t)
nagcurve <- lapply(nagcurve, melt)
nagcurve <- cbind(fit = nagcurve$fit[, -1], se = nagcurve$se.fit$value)
nagcurve[, 1] <- (100 - nagcurve[, 1]) / 100

ergcurve <- predict(ergen.wei[[1]], 
                    type = 'quantile',
                    p = seq(0.01, 0.99, by = 0.01),
                    se.fit = TRUE)
ergcurve <- lapply(ergcurve, function(x) x[1, ])
ergcurve <- lapply(ergcurve, t)
ergcurve <- lapply(ergcurve, melt)
ergcurve <- cbind(fit = ergcurve$fit[, -1], se = ergcurve$se.fit$value)
ergcurve[, 1] <- (100 - ergcurve[, 1]) / 100

# combine
reg.cur <- rbind(cbind(nacurve, loc = rep('NA', nrow(nacurve))),
                 cbind(ercurve, loc = rep('Eur', nrow(ercurve))))
greg.cur <- rbind(cbind(nagcurve, loc = rep('NA', nrow(nagcurve))),
                  cbind(ergcurve, loc = rep('Eur', nrow(ergcurve))))

regions <- rbind(cbind(reg.cur, heir = rep('species', nrow(reg.cur))),
                 cbind(greg.cur, heir = rep('genera', nrow(greg.cur))))

reg <- ggplot(regions, aes(x = fit.Var2, y = fit.value, fill = loc))
reg <- reg + geom_line(aes(colour = loc), size = 1)
reg <- reg + geom_ribbon(aes(ymin = fit.value - se, ymax = fit.value + se),
                         alpha = 0.3)
reg <- reg + coord_flip()
reg <- reg + facet_grid(. ~ heir)
reg <- reg + labs(y = 'Time', x = 'P(T > t)')
reg <- reg + scale_x_continuous(trans = log10_trans())
reg <- reg + scale_color_manual(values = cbp[-1],
                                name = 'Dietary\nCategory')
reg <- reg + scale_fill_manual(values = cbp[-1],
                               name = 'Dietary\nCategory')
reg <- reg + theme(axis.title.y = element_text(angle = 0),
                   axis.text = element_text(size = 20),
                   axis.title = element_text(size = 23),
                   legend.text = element_text(size = 17),
                   legend.title = element_text(size = 19),
                   strip.text = element_text(size = 20))
ggsave(filename = '../doc/figure/para_reg.png', plot = reg,
       width = 15, height = 10)


# diet
# species
nd <- predict(na.exp[[2]], newdata = data.frame(diet = c('carni',
                                                         'herb',
                                                         'insect',
                                                         'omni')),
              type = 'quantile',
              p = seq(0.01, 0.99, by = 0.01),
              se.fit = TRUE)
rownames(nd$fit) <- rownames(nd$se.fit) <- c('carni',
                                             'herb',
                                             'insect',
                                             'omni')
nd <- lapply(nd, t)
nd <- lapply(nd, melt)
nd <- cbind(fit = nd$fit, se = nd$se.fit$value)
nd[, 1] <- (100 - nd[, 1]) / 100

ed <- predict(er.wei[[2]], newdata = data.frame(diet = c('carni',
                                                         'herb',
                                                         'insect',
                                                         'omni')),
              type = 'quantile',
              p = seq(0.01, 0.99, by = 0.01),
              se.fit = TRUE)
rownames(ed$fit) <- rownames(ed$se.fit) <- c('carni',
                                             'herb',
                                             'insect',
                                             'omni')
ed <- lapply(ed, t)
ed <- lapply(ed, melt)
ed <- cbind(fit = ed$fit, se = ed$se.fit$value)
ed[, 1] <- (100 - ed[, 1]) / 100

# genera 
ndg <- predict(nagen.exp[[2]], newdata = data.frame(diet = c('carni',
                                                             'herb',
                                                             'insect',
                                                             'omni')),
               type = 'quantile',
               p = seq(0.01, 0.99, by = 0.01),
               se.fit = TRUE)
rownames(ndg$fit) <- rownames(ndg$se.fit) <- c('carni',
                                               'herb',
                                               'insect',
                                               'omni')
ndg <- lapply(ndg, t)
ndg <- lapply(ndg, melt)
ndg <- cbind(fit = ndg$fit, se = ndg$se.fit$value)
ndg[, 1] <- (100 - ndg[, 1]) / 100

edg <- predict(ergen.wei[[2]], newdata = data.frame(diet = c('carni',
                                                             'herb',
                                                             'insect',
                                                             'omni')),
               type = 'quantile',
               p = seq(0.01, 0.99, by = 0.01),
               se.fit = TRUE)
rownames(edg$fit) <- rownames(edg$se.fit) <- c('carni',
                                               'herb',
                                               'insect',
                                               'omni')
edg <- lapply(edg, t)
edg <- lapply(edg, melt)
edg <- cbind(fit = edg$fit, se = edg$se.fit$value)
edg[, 1] <- (100 - edg[, 1]) / 100

# combine
d.cur <- rbind(cbind(nd, loc = rep('NA', nrow(nd))),
               cbind(ed, loc = rep('Eur', nrow(ed))))
gd.cur <- rbind(cbind(ndg, loc = rep('NA', nrow(ndg))),
                cbind(edg, loc = rep('Eur', nrow(edg))))

diets <- rbind(cbind(d.cur, heir = rep('species', nrow(d.cur))),
               cbind(gd.cur, heir = rep('genera', nrow(gd.cur))))

die <- ggplot(diets, aes(x = fit.Var1, y = fit.value, fill = fit.Var2))
die <- die + geom_line(aes(colour = fit.Var2), size = 1)
die <- die + geom_ribbon(aes(ymin = fit.value - se, ymax = fit.value + se),
                         alpha = 0.3)
die <- die + coord_flip()
die <- die + facet_grid(loc ~ heir)
die <- die + scale_x_continuous(trans = log10_trans())
die <- die + labs(y = 'Time', x = 'P(T > t)')
die <- die + scale_color_manual(values = cbp[-1],
                                name = 'Dietary\nCategory')
die <- die + scale_fill_manual(values = cbp[-1],
                               name = 'Dietary\nCategory')
die <- die + theme(axis.title.y = element_text(angle = 0),
                   axis.text = element_text(size = 20),
                   axis.title = element_text(size = 23),
                   legend.text = element_text(size = 17),
                   legend.title = element_text(size = 19),
                   strip.text = element_text(size = 20))
ggsave(filename = '../doc/figure/para_diet.png', plot = die,
       width = 15, height = 10)

# locomotor category 
# species
nl <- predict(na.exp[[3]], newdata = data.frame(move = c('arboreal',
                                                         'ground dwelling',
                                                         'scansorial')),
              type = 'quantile',
              p = seq(0.01, 0.99, by = 0.01),
              se.fit = TRUE)
rownames(nl$fit) <- rownames(nl$se.fit) <- c('arboreal',
                                             'ground dwelling',
                                             'scansorial')
nl <- lapply(nl, t)
nl <- lapply(nl, melt)
nl <- cbind(fit = nl$fit, se = nl$se.fit$value)
nl[, 1] <- (100 - nl[, 1]) / 100

el <- predict(er.wei[[3]], newdata = data.frame(move = c('arboreal',
                                                         'ground dwelling',
                                                         'scansorial')),
              type = 'quantile',
              p = seq(0.01, 0.99, by = 0.01),
              se.fit = TRUE)
rownames(el$fit) <- rownames(el$se.fit) <- c('arboreal',
                                             'ground dwelling',
                                             'scansorial')
el <- lapply(el, t)
el <- lapply(el, melt)
el <- cbind(fit = el$fit, se = el$se.fit$value)
el[, 1] <- (100 - el[, 1]) / 100

# genera 
nlg <- predict(nagen.exp[[3]], newdata = data.frame(move = c('arboreal',
                                                             'ground dwelling',
                                                             'scansorial')),
               type = 'quantile',
               p = seq(0.01, 0.99, by = 0.01),
               se.fit = TRUE)
rownames(nlg$fit) <- rownames(nlg$se.fit) <- c('arboreal',
                                               'ground dwelling',
                                               'scansorial')
nlg <- lapply(nlg, t)
nlg <- lapply(nlg, melt)
nlg <- cbind(fit = nlg$fit, se = nlg$se.fit$value)
nlg[, 1] <- (100 - nlg[, 1]) / 100

elg <- predict(ergen.wei[[3]], newdata = data.frame(move = c('arboreal',
                                                             'ground dwelling',
                                                             'scansorial')),
               type = 'quantile',
               p = seq(0.01, 0.99, by = 0.01),
               se.fit = TRUE)
rownames(elg$fit) <- rownames(elg$se.fit) <- c('arboreal',
                                               'ground dwelling',
                                               'scansorial')
elg <- lapply(elg, t)
elg <- lapply(elg, melt)
elg <- cbind(fit = elg$fit, se = elg$se.fit$value)
elg[, 1] <- (100 - elg[, 1]) / 100


l.cur <- rbind(cbind(nl, loc = rep('NA', nrow(nl))),
               cbind(el, loc = rep('Eur', nrow(el))))
gl.cur <- rbind(cbind(nlg, loc = rep('NA', nrow(nlg))),
                cbind(elg, loc = rep('Eur', nrow(elg))))

moves <- rbind(cbind(l.cur, heir = rep('species', nrow(l.cur))),
               cbind(gl.cur, heir = rep('genera', nrow(gl.cur))))

loc <- ggplot(moves, aes(x = fit.Var1, y = fit.value, fill = fit.Var2))
loc <- loc + geom_line(aes(colour = fit.Var2), size = 1)
loc <- loc + geom_ribbon(aes(ymin = fit.value - se, ymax = fit.value + se),
                         alpha = 0.3)
loc <- loc + coord_flip()
loc <- loc + facet_grid(loc ~ heir)
loc <- loc + scale_x_continuous(trans = log10_trans())
loc <- loc + labs(y = 'Time', x = 'P(T > t)')
loc <- loc + scale_color_manual(values = cbp[-1],
                                name = 'Locomotor\nCategory')
loc <- loc + scale_fill_manual(values = cbp[-1],
                               name = 'Locomotor\nCategory')
loc <- loc + theme(axis.title.y = element_text(angle = 0),
                   axis.text = element_text(size = 20),
                   axis.title = element_text(size = 23),
                   legend.text = element_text(size = 17),
                   legend.title = element_text(size = 19),
                   strip.text = element_text(size = 20))
ggsave(filename = '../doc/figure/para_move.png', plot = loc,
       width = 15, height = 10)

# mass
# species
masses <- data.frame(mass = c(min(na.ecol$mass),
                              quantile(na.ecol$mass, .25),
                              quantile(na.ecol$mass, .5),
                              quantile(na.ecol$mass, .75),
                              max(na.ecol$mass)))
nm <- predict(na.exp[[4]], newdata = masses,
              type = 'quantile',
              p = seq(0.01, 0.99, by = 0.01),
              se.fit = TRUE)
rownames(nm$fit) <- rownames(nm$se.fit) <- c('Min', 'Q1', 'Median', 'Q3', 'Max')
nm <- lapply(nm, t)
nm <- lapply(nm, melt)
nm <- cbind(fit = nm$fit, se = nm$se.fit$value)
nm[, 1] <- (100 - nm[, 1]) / 100

# genera
gmasses <- data.frame(mass = c(min(na.genecol$mass),
                               quantile(na.genecol$mass, .25),
                               quantile(na.genecol$mass, .5),
                               quantile(na.genecol$mass, .75),
                               max(na.genecol$mass)))
nmg <- predict(nagen.exp[[4]], newdata = gmasses,
               type = 'quantile',
               p = seq(0.01, 0.99, by = 0.01),
               se.fit = TRUE)
rownames(nmg$fit) <- rownames(nmg$se.fit) <- c('Min', 'Q1', 'Median', 'Q3', 'Max')
nmg <- lapply(nmg, t)
nmg <- lapply(nmg, melt)
nmg <- cbind(fit = nmg$fit, se = nmg$se.fit$value)
nmg[, 1] <- (100 - nmg[, 1]) / 100

mas <- rbind(cbind(nm, heir = rep('species', nrow(nm))),
             cbind(nmg, heir = rep('genera', nrow(nmg))))
gma <- ggplot(mas, aes(x = fit.Var1, y = fit.value, fill = fit.Var2))
gma <- gma + geom_line(aes(colour = fit.Var2), size = 1)
gma <- gma + geom_ribbon(aes(ymin = fit.value - se, ymax = fit.value + se),
                         alpha = 0.3)
gma <- gma + coord_flip()
gma <- gma + facet_grid(. ~ heir)
gma <- gma + labs(y = 'Time', x = 'P(T > t)')
gma <- gma + scale_x_continuous(trans = log10_trans())
gma <- gma + scale_color_manual(values = cbp[-1],
                                name = 'Quantile\nlog(Mass)')
gma <- gma + scale_fill_manual(values = cbp[-1],
                               name = 'Quantile\nlog(Mass)')
gma <- gma + theme(axis.title.y = element_text(angle = 0),
                   axis.text = element_text(size = 20),
                   axis.title = element_text(size = 23),
                   legend.text = element_text(size = 17),
                   legend.title = element_text(size = 19),
                   strip.text = element_text(size = 20))
ggsave(filename = '../doc/figure/para_mass.png', plot = gma,
       width = 15, height = 10)


# climate
# species
nclim <- data.frame(climate = c(min(na.ecol$climate),
                                quantile(na.ecol$climate, .25),
                                quantile(na.ecol$climate, .5),
                                quantile(na.ecol$climate, .75),
                                max(na.ecol$climate)))
ncli <- predict(na.exp[[5]], newdata = nclim,
                type = 'quantile',
                p = seq(0.01, 0.99, by = 0.01),
                se.fit = TRUE)
rownames(ncli$fit) <- rownames(ncli$se.fit) <- c('Min', 'Q1', 'Median', 'Q3', 'Max')
ncli <- lapply(ncli, t)
ncli <- lapply(ncli, melt)
ncli <- cbind(fit = ncli$fit, se = ncli$se.fit$value)
ncli[, 1] <- (100 - ncli[, 1]) / 100

eclim <- data.frame(climate = c(min(er.ecol$climate),
                                quantile(er.ecol$climate, .25),
                                quantile(er.ecol$climate, .5),
                                quantile(er.ecol$climate, .75),
                                max(er.ecol$climate)))
ecli <- predict(er.wei[[4]], newdata = eclim,
                type = 'quantile',
                p = seq(0.01, 0.99, by = 0.01),
                se.fit = TRUE)
rownames(ecli$fit) <- rownames(ecli$se.fit) <- c('Min', 'Q1', 'Median', 'Q3', 'Max')
ecli <- lapply(ecli, t)
ecli <- lapply(ecli, melt)
ecli <- cbind(fit = ecli$fit, se = ecli$se.fit$value)
ecli[, 1] <- (100 - ecli[, 1]) / 100

# genera
ngclim <- data.frame(climate = c(min(na.genecol$climate),
                                 quantile(na.genecol$climate, .25),
                                 quantile(na.genecol$climate, .5),
                                 quantile(na.genecol$climate, .75),
                                 max(na.genecol$climate)))
ngcli <- predict(nagen.exp[[5]], newdata = ngclim,
                 type = 'quantile',
                 p = seq(0.01, 0.99, by = 0.01),
                 se.fit = TRUE)
rownames(ngcli$fit) <- rownames(ngcli$se.fit) <- c('Min', 'Q1', 'Median', 'Q3', 'Max')
ngcli <- lapply(ngcli, t)
ngcli <- lapply(ngcli, melt)
ngcli <- cbind(fit = ngcli$fit, se = ngcli$se.fit$value)
ngcli[, 1] <- (100 - ngcli[, 1]) / 100

egclim <- data.frame(climate = c(min(er.genecol$climate),
                                 quantile(er.genecol$climate, .25),
                                 quantile(er.genecol$climate, .5),
                                 quantile(er.genecol$climate, .75),
                                 max(er.genecol$climate)))
egcli <- predict(ergen.wei[[4]], newdata = egclim,
                 type = 'quantile',
                 p = seq(0.01, 0.99, by = 0.01),
                 se.fit = TRUE)
rownames(egcli$fit) <- rownames(egcli$se.fit) <- c('Min', 'Q1', 'Median', 'Q3', 'Max')
egcli <- lapply(egcli, t)
egcli <- lapply(egcli, melt)
egcli <- cbind(fit = egcli$fit, se = egcli$se.fit$value)
egcli[, 1] <- (100 - egcli[, 1]) / 100

# combine
c.cur <- rbind(cbind(ncli, loc = rep('NA', nrow(ncli))),
               cbind(ecli, loc = rep('Eur', nrow(ecli))))
gc.cur <- rbind(cbind(ngcli, loc = rep('NA', nrow(ngcli))),
                cbind(egcli, loc = rep('Eur', nrow(egcli))))

cli <- rbind(cbind(c.cur, heir = rep('species', nrow(c.cur))),
             cbind(gc.cur, heir = rep('genera', nrow(gc.cur))))


gcl <- ggplot(cli, aes(x = fit.Var1, y = fit.value, fill = fit.Var2))
gcl <- gcl + geom_line(aes(colour = fit.Var2), size = 1)
gcl <- gcl + geom_ribbon(aes(ymin = fit.value - se, ymax = fit.value + se),
                         alpha = 0.3)
gcl <- gcl + coord_flip()
gcl <- gcl + facet_grid(loc ~ heir)
gcl <- gcl + labs(y = 'Time', x = 'P(T > t)')
gcl <- gcl + scale_x_continuous(trans = log10_trans())
gcl <- gcl + scale_color_manual(values = cbp[-1],
                                name = 'Quantile\ndelta O^18')
gcl <- gcl + scale_fill_manual(values = cbp[-1],
                               name = 'Quantile\ndelta O^18')
gcl <- gcl + theme(axis.title.y = element_text(angle = 0),
                   axis.text = element_text(size = 20),
                   axis.title = element_text(size = 23),
                   legend.text = element_text(size = 17),
                   legend.title = element_text(size = 19),
                   strip.text = element_text(size = 20))
ggsave(filename = '../doc/figure/para_clim.png', plot = gcl,
       width = 15, height = 10)
