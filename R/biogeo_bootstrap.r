#' Resample biogeographic network statistic
#'
#' Process similar to that described in Sidor et al. 2013. 
#' This only performs a single resample.
#'
#' TODO add the bipartite check
#'
#' @param graph object of class igraph (bipartite)
#' @param fun function describing the biogeographic structure of the network
#' @param taxon hierarchical taxonomic information for each observation
#' @param data PBDB taxon list
#' @return
#' @export
#' @keywords
#' @author Peter D Smits <psmits@uchicago.edu>
#' @references
#' @examples
biogeo.re <- function(graph, fun, taxon, data) {
  bip <- bipartite.projection(graph)
  len <- lapply(bip, function(x) vcount(x))
  wt <- unlist(lapply(bip, function(x) all(V(x)$name %in% data$name.bi)))
  taxa <- V(bip[wt][[1]])$name
  loc <- V(bip[!wt][[1]])$name

  l.small <- ifelse(length(taxa) > length(loc), TRUE, FALSE)

  # remove a taxon at some probability
  # add back in a taxon of the same taxonomic group

  tot <- ecount(graph)
  # get neighbors of each locality
  pres <- abse <- list()
  for(ii in seq(length(loc))) {
    pres[[ii]] <- neighbors(graph, loc[ii])
    #abse[[ii]] <- seq(length(taxa))[-pres[[ii]]]
  }

  occ <- table(get.edgelist(graph)[, 1])
  prob.occ <- occ / sum(occ)

  rms <- lapply(pres, function(x, y) {
                oo <- y[x]
                nn <- ifelse(runif(length(oo)) < oo, TRUE, FALSE)
                x[nn]}, y = prob.occ)
  
  # taxonomic group info for rms
  tax <- lapply(rms, function(x, y) {
                na <- y[x]}, y = taxon)
  # grab a random member of the same taxonomic group
  spl.tax <- split(taxa, taxon)

  wsp <- lapply(tax, function(x, y) which(names(y) %in% x),
               y = spl.tax)

  ads <- vector(mode = 'list', length = length(tax))
  for(ii in seq(length(tax))) {
    if(length(tax[[ii]]) > 0) {
      ads[[ii]] <- laply(wsp[[ii]], function(x) {
                         sample(spl.tax[[x]], 1)})
    }
  }
  ads.n <- lapply(ads, function(x) which(V(graph)$name %in% x))

  # remove bad ids
  for(ii in seq(length(loc))) {
    ll <- which(V(graph)$name == loc[ii])
    rr <- cbind(rms[[ii]], ll)
    tr <- apply(rr, 1, function(x) E(graph, path = x))
    graph <- delete.edges(graph, tr)
  }

  # add new ids
  for(jj in seq(length(loc))) {
    ll <- which(V(graph)$name == loc[jj])
    uu <- cbind(ads.n[[jj]], ll)
    if(length(ads.n[[jj]] != 0)) graph <- add.edges(graph, as.vector(t(uu)))
  }

  # get rid of any taxa that have no neighbors
  pp <- lapply(taxa, function(x) neighbors(graph, x))
  rn <- which(V(graph)$name %in% taxa[which(unlist(lapply(pp, length)) == 0)])
  graph <- delete.vertices(graph, rn)

  # biogeographic summary statistic
  out <- fun(graph, l.small = l.small)

  out
}

#' Bootstrap biogeographic network statistic
#'
#' Process similar to that described in Sidor et al. 2013. 
#' This produces the full bootstrap distribution
#'
#' TODO add the bipartite check
#'
#' @param graph object of class igraph (bipartite)
#' @param fun function describing the biogeographic structure of the network
#' @param taxon hierarchical taxonomic information for each observation
#' @param data PBDB taxon list
#' @param nsim number of bootstrap replicates
#' @return
#' @export
#' @keywords
#' @author Peter D Smits <psmits@uchicago.edu>
#' @references
#' @examples
biogeo.boot <- function(graph, fun, taxon, data, nsim = 1000){
  out <- array(dim = nsim)
  for(ii in seq(nsim)) {
    out[ii] <- biogeo.re(graph = graph, 
                         fun = fun, 
                         taxon = taxon,
                         data = data)
  }

  out
}
