library(ape)
#source('http://www.phytools.org/read.newick/v0.4/read.newick.R')

#' convert a taxonomy into a phylogeny
#'
#' @param taxa df; 4 column data frame (genus, family, order, binomial)
#' @return object of class phylo
make.tree <- function(taxo) {
  taxo$occurrence.genus_name <- as.factor(taxo$occurrence.genus_name)
  taxo$family_name <- as.factor(taxo$family_name)
  taxo$order_name <- as.factor(taxo$order_name)
  taxo$name.bi <- as.factor(taxo$name.bi)
  tree <- as.phylo.formula(~ order_name/family_name/occurrence.genus_name/name.bi,
                           data = taxo)
  tree
}

#' make a tree of all the taxa
#'
#' @param data df; needs order, family, genus, binomial name
#' @return object of class phylo, all edges length 1
#' @export
big.tree <- function(data) {
  taxa <- unique(data[, c('order_name', 'family_name', 
                          'occurrence.genus_name', 'name.bi')])
  big.tree <- make.tree(taxa)
  big.tree$edge.length <- rep(1, nrow(big.tree$edge))
  big.tree
}

#' get the taxa for every locality
#'
#' @param graph bipartite graph object
#' @param l.small logical; is the smaller part of the bipartite the localities?
#' @return list of occurrences per site
#' @export
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
