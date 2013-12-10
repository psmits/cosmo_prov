library(ape)
#source('http://www.phytools.org/read.newick/v0.4/read.newick.R')

# convert a taxonomy into a phylogeny
make.tree <- function(taxo) {
  taxo$occurrence.genus_name <- as.factor(taxo$occurrence.genus_name)
  taxo$family_name <- as.factor(taxo$family_name)
  taxo$order_name <- as.factor(taxo$order_name)
  taxo$name.bi <- as.factor(taxo$name.bi)
  tree <- as.phylo.formula(~ order_name/family_name/occurrence.genus_name/name.bi,
                           data = taxo)
  tree
}

# make a tree of all the taxa
big.tree <- function(data) {
  taxa <- unique(data[, c('order_name', 'family_name', 
                          'occurrence.genus_name', 'name.bi')])
  big.tree <- make.tree(taxa)
  big.tree$edge.length <- rep(1, nrow(big.tree$edge))
  big.tree
}

# get the taxa for every locality
node.taxo <- function(graph, l.small = TRUE) {
  bip <- bipartite.projection(graph)
  len <- lapply(bip, function(x) length(V(x)))
  ws <- which.min(unlist(len))
  wm <- which.max(unlist(len))

  if(l.small) {
    st <- V(bip[[ws]])$name
    tx <- V(bip[[wm]])$name
  } else if(!l.small) {
    st <- V(bip[[wm]])$name
    tx <- V(bip[[ws]])$name
  }

  pres <- lapply(st, function(x) neighbors(graph, x))
  
  occ <- lapply(pres, function(x) tx[x])
  occ                                  
}

# grab the necessary tips, get mean cophenetic
mean.coph <- function(graph, data, l.small = TRUE) {
  tt <- node.taxo(graph, l.small)
  bt <- big.tree(data)

  co <- combn(length(tt), m = 2)
  ph <- array(dim = ncol(co))
  for(ii in seq(ncol(co))) {
    tip <- c(tt[[co[1, ii]]], tt[[co[2, ii]]])
    if(Reduce(identical, tip)) {
      ph[ii] <- 0
    } else {
      tr <- drop.tip(bt, which(!(bt$tip.label %in% tip)))
      cop <- cophenetic(tr)
     ph[ii] <- mean(cop[lower.tri(cop)])
    }
  }

  ph
}
