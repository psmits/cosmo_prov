#' Remove uncertainty from taxonomic naming
#'
#' This function takes an occurence matrix, as from the PBDB, with taxon names
#' that might be under some kind of uncertainty. This method then collapses 
#' the uncertainty.
#'
#' @param pbdb object of class pbdb
#' @param level regular expression identifying the exact column to collapse to
#'        fails of more than one result
#' @keywords
#' @export
#' @examples
collapse.names <- function(pbdb, level = 'genus_name') {
  if(!(class(pbdb) == 'pbdb')) {
    stop('object is not of class pbdb')
  }

  occ <- pbdb$occurence

  w.sites <- grep(pattern = 'X[0-9]', x = colnames(occ), perl = TRUE)
  w.taxa <- grep(pattern = level, x = colnames(occ), perl = TRUE)
  taxa <- occ[, w.taxa]
  sites <- occ[, w.sites]

  # if level is at 'species', make sure to preserve bionmial
  # need to also collapse the ecological information
  occ <- collapse.taxa(taxa, sites, name = colnames(occ)[w.taxa])
  
  pbdb$occurence <- occ
 
  if(is.element('ecology', names(pbdb))) {
    # clean ecologies
    dup <- duplicated(taxa)
    pbdb$ecology <- pbdb$ecology[!dup, ]
  }
  if(is.element('geology', names(pbdb))) {
    pbdb$geology <- pbdb$geology[, c(w.taxa, w.sites)]
  }

  pbdb
}

#' Collapse same named things
#'
#' Internal function to support rm.uncer
#'
#' TODO problems with mismatch lengths that don't make sense
#'
#' @param taxa name information from a pbdb$occurence
#' @param sites abundance information
collapse.taxa <- function(taxa, sites, name) {
  require(plyr)
  ab <- ddply(cbind(taxa, sites), .(taxa), function(x) colSums(x[, -1]))

  colnames(ab)[1] <- name
  ab
}
