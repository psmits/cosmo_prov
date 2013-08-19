#' Clean a Paleobiological Database occurence matrix
#'
#' This function takes the raw output of a Paleobiological Database occurence 
#' matrix and seperates additional information that exists outside of the specific
#' multivariate abundance matrix. This includes ecological information for the 
#' taxa and geological information for the sites.
#'
#' TODO include logicals if there is, or isn't, ecological or geological
#' information
#'
#' @param raw.pbdb raw output occurence matrix from the PBDB
#' @keywords
#' @export
#' @examples
clean.pbdb <- function(raw.pbdb) {
  out <- list()
  cc <- colnames(raw.pbdb)

  # i'm lucky that sites begin with an X and are followed by numbers
  sites <- grep(pattern = 'X[0-9]', x = cc, perl = TRUE)

  # include logicals if there isn't collections information
  fst <- mam[, 1]
  geo <- grep(pattern = 'collections', x = fst, perl = TRUE)

  # include logicals if there isn't ecological information
  # occurence information is between the top and the geology
  occ <- raw.pbdb[1:(min(geo) - 1), 1:(max(sites))]
  # ecology information is after the sites
  eco <- raw.pbdb[-geo, max(sites):ncol(mam)]
  # geological information is at the bottom of the occurence information
  collec <- raw.pbdb[geo, 1:(max(sites))]

  # break genus from occurence
  # make occurence numeric, then put back together
  oo <- occ[, sites]
  oo <- suppressWarnings(apply(oo, 2, function(x) as.numeric(as.character(x))))
  oo[is.na(oo)] <- 0
  occ[, sites] <- oo

  out$occurence <- occ
  out$ecology <- eco
  out$geology <- collec

  class(out) <- 'pbdb'
  out
}
