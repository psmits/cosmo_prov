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


#' Grab heir
#'
#'
grab.heir <- function(taxon, key) {
  eol.id <- my.get_eolid(taxon, key = key)
  tax <- classification(eol.id, key = key)
  names(tax) <- taxon

  oandf <- lapply(tax, function(x) {
                  if(class(x) != 'character') {
                    pick <- which(x$taxonRank == 'order')
                    x[seq(pick, nrow(x)), ]}})
  oandf <- oandf[!laply(oandf, is.null)]

  oandf <- lapply(oandf, function(x) {
                  data.frame(lapply(x, as.character), 
                             stringsAsFactors = FALSE)})
  o <- list()
  for(ii in seq(length(oandf))) {
    o[[ii]] <- rbind(oandf[[ii]][, 2:3], c(names(oandf)[ii], 'genus'))
  }
  names(o) <- names(oandf)

  o <- lapply(o, function(x) {
              x <- t(x)
              colnames(x) <- x[2, ]
              x <- x[1, ]
              x})
  o <- Reduce(rbind, o)
  rownames(o) <- NULL
  o
}
