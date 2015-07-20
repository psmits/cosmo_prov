library(ape)
library(plyr)
library(phytools)
library(geiger)
library(stringr)
library(paleotree)
library(phangorn)

load('../data/setup_tree.rdata')

#big.tree <- mrp.supertree(species.trees)
#save(big.tree, file = '../data/super_tree.rdata')
load('../data/super_tree.rdata')

# get rid of the stupid tips
if(class(big.tree) == 'multiPhylo') {
  spt <- big.tree[[1]] 
} else {
  spt <- big.tree
}

fad <- ddply(dat, .(name.bi), summarize, max(bins))
lad <- ddply(dat, .(name.bi), summarize, min(bins))
fad[, 1] <- str_replace(fad[, 1], ' ', '_')
lad[, 1] <- str_replace(lad[, 1], ' ', '_')
datmat <- cbind(fad[, 2], lad[, 2])
rownames(datmat) <- fad[, 1]

check <- name.check(spt, datmat)
spt <- drop.tip(spt, check$tree_not_data)
datmat <- datmat[!(rownames(datmat) %in% check$data_not_tree), ]
datmat <- datmat[match(spt$tip.label, rownames(datmat)), ]

spt <- timeLadderTree(spt, timeData = datmat)
spt <- timePaleoPhy(spt, timeData = datmat, 
                         type = 'mbl', vartime = 0.1)
save(spt, file = '../data/scaled_super.rdata')
