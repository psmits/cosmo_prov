#' Bionomial name maker
#'
#' Convenience function to make a bionomial name from seperate
#' genus and species information.
#'
#' @param genus string genus name
#' @param spec string species name
#' @param sep seperator between genus and species text
#' @return
#' @export
#' @author Peter D Smits <psmits@uchicago.edu>
#' @examples
binom.make <- function(genus, species, sep = ' ') {
  paste(genus, species, sep = sep)
}


#' Hierarchy finder
#'
#' Get the taxonomic hierarchy information for a taxa given binomial name.
#' This is an internal function that really only works with a particular 
#' data structure used in this study. Incidentally, that data structure is 
#' just a cleaned version of the PBDB taxon list.
#'
#' @param name binomial name such as the ouput from binom.make
#' @param level string column name in data with appropriate hierarchical information
#' @param data PBDB taxon list
#' @return
#' @export
#' @author Peter D Smits <psmits@uchicago.edu>
#' @examples
hierfind <- function(name, level, data) {
  nec <- data[, level]
  unique(nec[data$name.bi %in% name])
}


#' Get taxon info from network
#'
#' From a given bipartite network, recover taxonomic info for each taxon.
#'
#' @param graph object of class igraph (bipartite)
#' @param level string column name in data with appropriate hierarchical information
#' @param data PBDB taxon list
#' @return
#' @export
#' @author Peter D Smits <psmits@uchicago.edu>
#' @examples
get.hier <- function(graph, level, data) {
  bip <- bipartite.projection(graph)
  len <- lapply(bip, function(x) vcount(x))
  wt <- unlist(lapply(bip, function(x) all(V(x)$name %in% data$name.bi)))
  taxa <- V(bip[wt][[1]])$name
  loc <- V(bip[!wt][[1]])$name

  laply(taxa, hierfind, level = level, data = data)
}


#' replace taxonomy for data
#'
#' @param data; PBDB output
#' @param update; new hierarchy information (order, family, genus)
#' @return updated data object
#' @export
replace.taxonomy <- function(data, update) {
  for(ii in seq(nrow(update))) {
    rp <- which(data$occurrence.genus_name == update[ii, 3])
    data[rp, 'order_name'] <- update[ii, 1]
    data[rp, 'family_name'] <- update[ii, 2]
    rfa <- which(data$family_name == update[ii, 2])
    data[rfa, 'order_name'] <- update[ii, 1]
  }
  data
}


#' double check missing orders
#'
#' @param data; PBDB output
#' @return vector character family names
#' @export
miss.ord <- function(data) {
  mat.fam <- ord.fam(data)
  fix.ord <- tax_name(query = mat.fam, get = 'order')
  ups <- mat.fam[!is.na(fix.ord)]
  for(ii in seq(length(ups))) {
    rr <- data$family_name == ups[ii]
    data$order_name[rr] <- na.omit(fix.ord)[ii, 1]
  }
  data
}

#' families of missing orders
#'
#' @param data; PBDB output
#' @return vector of character family names
#' @export
ord.fam <- function(data) {
  ord <- data$order_name == ''
  mat.fam <- unique(data$family_name[ord])
  mat.fam <- mat.fam[mat.fam != '']
  mat.fam
}
