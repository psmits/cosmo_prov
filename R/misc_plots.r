library(ggplot2)
library(reshape2)
library(scales)

source('../R/surv_selection.r')

theme_set(theme_bw())
cbp <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
         "#D55E00", "#CC79A7")

# occupancy
spec.occ <- rbind(cbind(na.taxa.occ, loc = rep('NA', nrow(na.taxa.occ))),
                  cbind(er.taxa.occ, loc = rep('Eur', nrow(er.taxa.occ))))
gen.occ <- rbind(cbind(nag.taxa.occ, loc = rep('NA', nrow(nag.taxa.occ))),
                 cbind(erg.taxa.occ, loc = rep('Eur', nrow(erg.taxa.occ))))

occup <- rbind(cbind(spec.occ, heir = rep('species', nrow(spec.occ))),
               cbind(gen.occ, heir = rep('genera', nrow(gen.occ))))

go <- ggplot(occup, aes(x = occ))
go <- go + geom_histogram(binwidth = 1)
go <- go + labs(x = 'Mean BU occupancy')
go <- go + facet_grid(loc ~ heir)
go <- go + theme(axis.title.y = element_text(angle = 0),
                 axis.text = element_text(size = 20),
                 axis.title = element_text(size = 23),
                 legend.text = element_text(size = 17),
                 legend.title = element_text(size = 19),
                 strip.text = element_text(size = 20))
ggsave(filename = '../doc/figure/occ_dist.png', plot = go,
       width = 15, height = 10)
