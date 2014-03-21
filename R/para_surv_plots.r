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
