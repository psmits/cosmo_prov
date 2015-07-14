library(ape)
library(plyr)
library(phytools)
library(stringr)
library(paleotree)
library(phangorn)

load('../data/setup_tree.rdata')

big.tree <- mrp.supertree(species.trees)

fad <- ddply(dat, .(name.bi), summarize, max(bins))
lad <- ddply(dat, .(name.bi), summarize, min(bins))
fad[, 1] <- str_replace(fad[, 1], ' ', '_')
lad[, 1] <- str_replace(lad[, 1], ' ', '_')
fad <- fad[match(big.tree$tip.label, fad[, 1]), ]
lad <- lad[match(big.tree$tip.label, lad[, 1]), ]
datmat <- cbind(fad[, 2], lad[, 2])
rownames(datmat) <- fad[, 1]

# get rid of the stupid tips
if(class(big.tree) == 'multiPhylo') {
  spt <- big.tree[[1]] 
} else {
  spt <- big.tree
}

dr <- spt$tip.label[!(spt$tip.label %in% rownames(datmat))]
spt <- drop.tip(spt, dr)
spt <- timeLadderTree(spt, timeData = datmat)
spt <- timePaleoPhy(spt, timeData = datmat, 
                         type = 'mbl', vartime = 0.1)
save(spt, file = '../data/scaled_super.rdata')
