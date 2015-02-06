library(plyr)
library(phytools)
library(stringr)
library(paleotree)
library(phangorn)

load('../data/update_taxonomy.rdata')
source('../R/phylo_gen.r')
source('../R/na_mung.r')
#load('../data/many_trees.rdata')

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
new.tax$order_name[new.tax$order_name == ''] <- 'Unk'
new.tax$family_name[new.tax$family_name == ''] <- 'Unk'
for(ii in seq(nrow(new.tax))) {
  if(new.tax$order_name[ii] == 'Artiodactyla') {
    new.tax$order_name[ii] <- 'Cetartiodactyla'
  }
  if(new.tax$family_name[ii] == 'Unk') {
    new.tax$family_name[ii] <- paste0(new.tax$order_name[ii], 
                                      new.tax$family_name[ii])
  }
}

my.taxonomy <- new.tax[, c('order_name', 'family_name', 
                           'occurrence.genus_name', 'name.bi')]
my.taxonomy <- unique(my.taxonomy)
na.tree <- make.tree(my.taxonomy)
na.tree$tip.label <- str_replace(na.tree$tip.label, ' ', '_')
fad <- ddply(dat, .(name.bi), summarize, max(bins))
lad <- ddply(dat, .(name.bi), summarize, min(bins))
fad[, 1] <- str_replace(fad[, 1], ' ', '_')
lad[, 1] <- str_replace(lad[, 1], ' ', '_')
fad <- fad[match(na.tree$tip.label, fad[, 1]), ]
lad <- lad[match(na.tree$tip.label, lad[, 1]), ]
datmat <- cbind(fad[, 2], lad[, 2])
rownames(datmat) <- fad[, 1]

na.tree <- timeLadderTree(na.tree, timeData = datmat)
na.scale <- timePaleoPhy(na.tree, timeData = datmat, 
                         type = 'mbl', vartime = 0.1)
save(na.scale, file = '../data/taxonomy_tree.rdata')

# raia.tree is species level
# super tree is species level
species.trees <- c(raia.tree, super.tree[[1]], na.tree)
class(species.trees) <- 'multiPhylo'
big.tree <- mrp.supertree(species.trees)
save(species.trees, file = '../data/super_big_tree.rdata')

# get rid of the stupid tips
if(class(species.trees) == 'multiPhylo') {
  spt <- species.trees[[1]] 
} else {
  spt <- species.trees
}

dr <- spt$tip.label[!(spt$tip.label %in% dat$name.bi)]
spt <- drop.tip(spt, dr)
spt <- timeLadderTree(spt, timeData = datmat)
spt <- timePaleoPhy(spt, timeData = datmat, 
                         type = 'mbl', vartime = 0.1)
save(spt, file = '../data/scaled_super.rdata')



#to.genera <- function(tree) {
#  #  from liam revell
#  tips<-tree$tip.label
#  genera<-unique(sapply(strsplit(tips,"_"),function(x) x[1]))
#  ## here are our genera
#  ## now drop all but one of each
#  ii<-sapply(genera,function(x,y) grep(x,y)[1],y=tips)
#  tree<-drop.tip(tree,setdiff(tree$tip.label,tips[ii]))
#  tree$tip.label<-sapply(strsplit(tree$tip.label,"_"),function(x) x[1])
#  tree
#}
#
#not.clean <- laply(lapply(mam.flat, function(x) x$tip.label), is.bad)
#to.clean <- mam.flat[not.clean]
#
#scrub <- function(tree) {
#  tree <- collapse.singles(tree)
#  if(any(is.na(tree$tip.label)) | any(tree$tip.label == '')) return(NULL)
#  tester <- drop.tip(tree, which(str_detect(tree$tip.label, '\\.')))
#  if(sum(str_detect(tester$tip.label, '_')) == 1) {
#    return(NULL)
#  } else if (is.null(tester)) { 
#    return(NULL)
#  } else {
#    tester <- drop.tip(tester, which(!str_detect(tester$tip.label, '_')))
#  }
#  if(!is.null(tester)) {
#     tester$tip.label <- str_replace(tester$tip.label, '^([^_]*_[^_]*)_.*$', '\\1')
#  }
#  tester
#}
#
#clean.out <- llply(to.clean, scrub)
#class(mam.flat) <- 'NULL'
#mam.flat[not.clean] <- clean.out
#rmd <- which(laply(mam.flat, is.null))
#final.trees <- mam.flat[-rmd]
#no.extinct <- laply(final.trees, function(x) {
#                    nar <- clean.taxon(x$tip.label)
#                    o <- nar %in% unique(dat$name.bi)
#                    any(o)})
#
#my.trees <- final.trees[no.extinct]
#my.trees[[length(my.trees) + 1]] <- raia.tree
#
#
#my.trees[[length(my.trees) + 1]] <- na.tree
#
#genera.trees <- list()
#for(ii in seq(length(my.trees))) {
#  genera.trees[[ii]] <- try(to.genera(my.trees[[ii]]))
#}
#genera.trees <- genera.trees[laply(genera.trees, 
#                                   function(x) class(x) != 'try-error')]
#genera.trees[[length(genera.trees) + 1]] <- tom.tree
#
#class(my.trees) <- 'multiPhylo'
#class(genera.trees) <- 'multiPhylo'
##my.super <- mrp.supertree(my.trees)
##genera.super <- mrp.supertree(genera.trees)
