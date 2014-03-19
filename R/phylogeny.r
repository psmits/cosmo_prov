library(paleotree)
library(plyr)

source('../R/phylo_gen.r')

load('../data/update_taxonomy.rdata')

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

save(na.tree, er.tree,
     file = '../data/taxon_trees.rdata')
