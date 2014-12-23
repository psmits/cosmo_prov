library(igraph)

set.seed(420)
pdf(file = '../doc/metacommunity/figure/rng_graphs.pdf')
par(mfrow = c(3, 4), mar = c(1, 2, 2, 2))
for(i in seq(12)) plot(erdos.renyi.game(10, 0.2), 
                       vertex.label = NA, 
                       vertex.size = 30)
dev.off()
