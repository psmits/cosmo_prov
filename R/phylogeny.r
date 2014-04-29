library(paleotree)
library(plyr)
library(geiger)

source('../R/phylo_gen.r')

load('../data/update_taxonomy.rdata')

source('../R/na_mung.r')

# north america
no.fam <- which(na.tax$family_name == '')
no.ord <- which(na.tax$order_name == '')
rms <- unique(c(no.fam, no.ord))
clean.na <- na.tax[-rms, ]
na.tree <- big.tree(clean.na[, c('occurrence.genus_name', 
                                 'family_name', 
                                 'order_name', 
                                 'name.bi')])
time.data <- clean.na[, c('name.bi', 'ma_max', 'ma_min')]
time.data <- ddply(time.data, .(name.bi), summarize,
                   FAD = max(ma_max),
                   LAD = min(ma_min))
rownames(time.data) <- time.data[, 1]
time.data <- time.data[, 2:3]
na.tree <- timePaleoPhy(na.tree, time.data,
                        type = 'mbl', vartime = 0.01)


# europe
no.fam <- which(er.tax$family_name == '')
no.ord <- which(er.tax$order_name == '')
rms <- unique(c(no.fam, no.ord))
clean.er <- er.tax[-rms, ]
er.tree <- big.tree(clean.er[, c('occurrence.genus_name', 
                                 'family_name', 
                                 'order_name', 
                                 'name.bi')])

time.data <- clean.er[, c('name.bi', 'ma_max', 'ma_min')]
time.data <- ddply(time.data, .(name.bi), summarize,
                   FAD = max(ma_max),
                   LAD = min(ma_min))
rownames(time.data) <- time.data[, 1]
time.data <- time.data[, 2:3]
time.data[time.data[, 1] < time.data[, 2], ] <- 
  time.data[time.data[, 1] < time.data[, 2], 2:1]
er.tree <- timePaleoPhy(er.tree, time.data,
                        type = 'mbl', vartime = 0.01)


# make some plots
uti <- dat[, c('name.bi', 'comdiet', 'comlife')]
uti <- uti[!(duplicated(uti$name.bi)), ]
miss.match <- name.check(na.tree, data.names = uti$name.bi)
prune <- ape::drop.tip(na.tree, miss.match$tree_not_data)
sm <- uti[!(uti$name.bi %in% miss.match$data_not_tree), ]
sm <- sm[match(prune$tip.label, sm$name.bi), ]
sm$name.bi <- as.character(sm$name.bi)

dit <- sm$comdiet
cbp <- c('#E69F00', '#56B4E9', '#009E73', '#F0E442', 
         '#0072B2', '#D55E00', '#CC79A7')
dit[dit == 'carni'] <- cbp[1]
dit[dit == 'herb'] <- cbp[2]
dit[dit == 'insect'] <- cbp[3]
dit[dit == 'omni'] <- cbp[4]

pdf(file = '../doc/figure/na_phylo_diet.pdf')
par(mar = c(0, 0, 0, 0))
plot.phylo(prune, type = 'fan', show.tip.label = FALSE)
tiplabels(pch = 21, cex = 0.5, col = dit, bg = dit)
legend('bottomright', legend = unique(sm$comdiet),
       col = cbp[1:4], pch = 19, ncol = 1, cex = 1)
dev.off()

save(na.tree, er.tree,
     file = '../data/taxon_trees.rdata')
