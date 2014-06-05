library(ggplot2)
library(reshape2)

source('../R/surv_selection.r')

theme_set(theme_bw())
cbp <- c('#A6CEE3', '#B2DF8A', '#FB9a99', '#FDBF6F', '#CAB2D6', '#FFFF99',
         '#1F78B4', '#33A02C', '#E31A1C', '#FF7F00', '#6A3D9A', '#B15928')

imp <- rbind(cbind(na.imp, loc = rep('NA', nrow(na.imp))), 
             cbind(er.imp, loc = rep('Eur', nrow(er.imp))))
gimp <- ggplot(imp, aes(x = pred, y = imp, fill = loc))
gimp <- gimp + geom_bar(stat = 'identity', position = 'dodge')
gimp <- gimp + scale_fill_manual(values = cbp, name = 'Region')
gimp <- gimp + labs(x = 'variable', y = 'rel.\nimp.')
gimp <- gimp + theme(axis.title.y = element_text(angle = 0),
                     axis.text = element_text(size = 20),
                     axis.title = element_text(size = 23),
                     legend.text = element_text(size = 17),
                     legend.title = element_text(size = 19))
ggsave(filename = '../doc/figure/rel_imp.png', plot = gimp, 
       width = 15, height = 10)

