# miscellaneous functions for analyzing networks necessary

#' Split graph into list community subgraphs
#'
#' @param graph igraph object
#' @param vv vector of vertex memberships
#' @export
#' @keywords
#' @examples
list.subgraph <- function(graph, vv) {
  out <- list()
  for(ii in seq(max(vv$membership))) {
    out[[ii]] <- induced.subgraph(graph, which(vv$membership == ii))
  }
  out
}
