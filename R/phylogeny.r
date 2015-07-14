library(ape)
library(plyr)
library(phytools)
library(stringr)
library(paleotree)
library(phangorn)

load('../data/update_taxonomy.rdata')
source('../R/phylo_gen.r')
source('../R/na_mung.r')

raia.tree <- read.tree('../data/raia_tree.txt')
tom.tree <- read.nexus('../data/tomiya_tree.nex')

fam.tree <- read.nexus('../data/meredith.nex')
super.tree <- read.nexus('../data/bininda_emonds.nex')

# make taxonomy tree of what i have
# north america
no.fam <- which(na.tax$family_name == '')
no.ord <- which(na.tax$order_name == '')
rms <- unique(c(no.fam, no.ord))
clean.na <- na.tax[-rms, ]
uni.tax <- unique(clean.na[, 1:4])

new.tax <- replace.taxonomy(dat, uni.tax[, 1:3])
# some things don't have known orders and families. 
# so i want them to be polytomy at the "root"
new.tax$order_name[new.tax$order_name == 'Artiodactyla'] <- 'Cetartiodactyla'

# split by orders
by.order <- split(new.tax, new.tax$order_name)
order.tree <- list()
for(ii in seq(length(by.order))) {
  check <- by.order[[ii]]
  check <- check[, c('family_name', 'occurrence.genus_name', 'name.bi')]
  check$family_name <- as.factor(check$family_name)
  check$occurrence.genus_name <- as.factor(check$occurrence.genus_name)
  check$name.bi <- str_replace(check$name.bi, ' ', '_')
  check$name.bi <- as.factor(check$name.bi)
  check <- unique(check)

  # screen situations with only 1 family, causes problems
  if(length(levels(check[, 1])) == 1) {
    temp <- as.phylo.formula(~ occurrence.genus_name/
                             name.bi, 
                             data = check)
  } else {
    temp <- as.phylo.formula(~ family_name/
                             occurrence.genus_name/
                             name.bi, 
                             data = check)
  }
  temp <- collapse.singles(temp)

  # make the unknown family node "disappear"
  temp <- unitLengthTree(temp)
  roo <- getMRCA(temp, check$name.bi[check$family_name == ''])
  temp$edge.length[apply(temp$edge, 1, function(x) any(x == roo))] <- 0
  temp <- di2multi(temp)

  order.tree[[ii]] <- temp
}
names(order.tree) <- names(by.order)
ntip <- laply(order.tree, function(x) length(x$tip.label))
order.tree <- order.tree[ntip != 1]
# combine the orders into a huge tree
taxon.tree <- Reduce(bind.tree, order.tree)

species.trees <- list(raia.tree, super.tree[[1]], taxon.tree)
class(species.trees) <- 'multiPhylo'

save.image(file = '../data/setup_tree.rdata')
