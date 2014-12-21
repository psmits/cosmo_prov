library(rstan)
library(arm)
library(parallel)
library(xtable)
library(stringr)
library(plyr)
library(igraph)
library(ape)

RNGkind(kind = "L'Ecuyer-CMRG")
seed <- 420
nsim <- 1000

set.seed(seed)
source('../R/surv_setup.r')
load('../data/taxonomy_tree.rdata')

species.graph <- function(bipartite) {
  bip <- bipartite.projection(bipartite)
  if(any(grepl('[0-9]', V(bip[[1]])$name))) {
    taxa <- bip[[2]]
    st <- bip[[1]] 
  } else {
    taxa <- bip[[1]]
    st <- bip[[2]] 
  }
  out <- taxa
  out
}

spg <- llply(taxawin, function(x) {
             o <- species.graph(x)
             E(o)$weight <- 1
             o})
name <- llply(spg, function(x) V(x)$name)
inin <- llply(name, function(x) x %in% na.ecol$taxa)
spg <- Map(function(x, y) delete.vertices(x, which(!y)), spg, inin)

adj <- llply(spg, function(x) get.adjacency(x, sparse = FALSE))
deg <- llply(spg, degree)
size <- llply(spg, vcount)
name <- llply(spg, function(x) V(x)$name)

# have to make these line up correctly
ecol <- llply(name, function(x) {
              tt <- na.ecol[na.ecol$taxa %in% x, ]
              tt[match(x, tt[, 1]), ]
              tt})


# prepare the phylo vcv
tax.nam <- str_replace(na.ecol[, 1], ' ', '_')
to.drop <- na.scale$tip.label[!(na.scale$tip.label %in% tax.nam)]
na.tree <- drop.tip(na.scale, to.drop)

# for each window
vcv.s <- list()
for(ii in seq(length(name))) {
  nono <- na.tree$tip.label[!(na.tree$tip.label %in% 
                              str_replace(name[[ii]], ' ', '_'))]
  temp <- drop.tip(na.tree, nono)
  temp$edge.length <- temp$edge.length / max(diag(vcv(temp)))
  temp.vcv <- vcv(temp)
  
  # have to have this line up correctly
  o <- match(str_replace(name[[ii]], ' ', '_'), colnames(temp.vcv))
  temp.vcv <- temp.vcv[o, o]

  vcv.s[[ii]] <- temp.vcv  
}

