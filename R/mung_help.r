#' Assign lat long grid
#'
#' @param lat
#' @param long
#' @param width
#' @return
#' @export
#' @author Peter D Smits <psmits@uchicago.edu>
#' @references
#' @examples
grid.id <- function(lat, long, width = 2) {
  mlat <- max(lat)
  mlng <- min(long)

  plat <- seq(from = floor(min(lat)), to = ceiling(mlat) + 1, by = width)
  plng <- seq(from = ceiling(max(long)), to = floor(mlng) - 1, by = -width)

  glat <- cut(lat, breaks = plat, include.lowest = TRUE)
  glng <- cut(long, breaks = plng, include.lowest = TRUE)

  gid <- interaction(glat, glng)
  gid
}
