#' Good's estimate of frequency coverage
#'
#' Coverage is an estimator of the amount of observed diversity of some 
#' set of categorical variables. For example, given a distribution of 
#' taxonomic abundances, it is possible to determine how much of the possible 
#' set has been sampled.
#'
#' @param ab table of total observed abundances
#' @return
#' @export
#' @keywords
#' @author
#' @references
#' @examples
coverage <- function(ab) {
  oo <- sum(ab)
  ss <- sum(ab == 1)
  if (ss == oo) ss = oo - 1
  uu <- 1 - ss / oo
  uu
}
