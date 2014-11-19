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


extract.genus <- function(string) {
  g <- str_extract(string, '^[^_ ]*')
  g
}
clean.taxon <- function(string) {
  g <- str_replace(string, '[_ ]', ' ')
  g
}
ext.clade <- function(tree, label) {
  # extract a clade that actually works
  root <- length(tree$tip.label) + which(tree$node.label == label)
  clade <- getDescendants(tree, root)
  global.root <- length(tree$tip.label) + 1
  tips <- clade[clade < global.root]
  excl <- ape::drop.tip(tree, seq(global.root - 1)[-tips])
  excl
}
labeled.subclades <- function(tree, name, fam) {
  # grab a labeled sublcade of your choice
  nono <- which(!(clean.taxon(tree$tip.label) %in% name))
  small <- drop.tip(tree, nono)
  
  # node matches
  ord <- small$node.label[small$node.label != '']

  matchs <- sort(unique(c(fam[which(fam %in% ord)])))
  subs <- llply(matchs, function(x) ext.clade(small, x))
  names(subs) <- matchs
  subs
}
is.bad <- function(string) {
  # find the problem trees
  first <- any(str_detect(string, '[\\.0-9]'))

  if(first) {
    return(first) 
  } else {
    second <- any(!str_detect(string, '_'))
    return(second)
  }
}

