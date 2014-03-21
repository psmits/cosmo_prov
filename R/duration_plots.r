library(survival)
library(ggplot2)
library(scales)
library(reshape2)

source('../R/na_surv.r')
source('../R/eur_surv.r')

# north america
nadur$diet <- na.ecol$diet
nadur$move <- na.ecol$move
nadur <- nadur[order(nadur$fad), ]
nadur$num <- seq(nrow(nadur))
na.seg <- ggplot(nadur, aes(x = fad, y = num))
na.seg <- na.seg + geom_segment(aes(xend = lad, yend = num))
ggsave(filename = '../doc/figure/na_dur_seg.png', plot = na.seg,
       width = 15, height = 10)


erdur$diet <- er.ecol$diet
erdur$move <- er.ecol$move
erdur <- erdur[order(erdur$fad), ]
erdur$num <- seq(nrow(erdur))
er.seg <- ggplot(erdur, aes(x = fad, y = num))
er.seg <- er.seg + geom_segment(aes(xend = lad, yend = num))
ggsave(filename = '../doc/figure/er_dur_seg.png', plot = er.seg,
       width = 15, height = 10)
