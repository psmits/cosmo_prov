library(plyr)
library(phytools)
library(stringr)

load('../data/update_taxonomy.rdata')
source('../R/phylo_gen.r')
source('../R/na_mung.r')
load('../data/many_trees.rdata')

raia.tree <- read.tree('../data/raia_tree.txt')
tom.tree <- read.nexus('../data/tomiya_tree.nex')

not.clean <- laply(lapply(mam.flat, function(x) x$tip.label), is.bad)

to.clean <- mam.flat[not.clean]

scrub <- function(tree) {
  tree <- collapse.singles(tree)

  if(any(is.na(tree$tip.label)) | any(tree$tip.label == '')) return(NULL)

  tester <- drop.tip(tree, which(str_detect(tree$tip.label, '\\.')))
  if(sum(str_detect(tester$tip.label, '_')) == 1) {
    return(NULL)
  } else if (is.null(tester)) { 
    return(NULL)
  } else {
    tester <- drop.tip(tester, which(!str_detect(tester$tip.label, '_')))
  }

  if(!is.null(tester)) {
     tester$tip.label <- str_replace(tester$tip.label, '^([^_]*_[^_]*)_.*$', '\\1')
  }
  tester
}

clean.out <- llply(to.clean, scrub)

class(mam.flat) <- 'NULL'
mam.flat[not.clean] <- clean.out
rmd <- which(laply(mam.flat, is.null))

final.trees <- mam.flat[-rmd]
no.extinct <- laply(final.trees, function(x) {
                    nar <- clean.taxon(x$tip.label)
                    o <- nar %in% unique(dat$name.bi)
                    any(o)})

my.trees <- final.trees[no.extinct]
my.trees[[length(my.trees)]] <- raia.tree

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

my.trees[[length(my.trees)]] <- na.tree

class(my.trees) <- 'multiPhylo'
#my.super <- mrp.supertree(my.trees)
