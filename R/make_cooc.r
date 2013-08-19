require(reshape2)

#' Make co-occurrence from multivariate abundance.
#'
#' This function creates a co-occurrence matrix from a multivariate abundance
#' matrix. By default, this function assumes that rows are sites and columns
#' are taxa. This is end-user modifiable TODO.
#'
#' TODO testing
#'
#' @param mat multivariate abundance matrix
#' @param weighted is the resulting matrix supposed to be weighted by abundance?
#' @keywords
#' @export
#' @examples
make.cooc <- function(mat, weighted = FALSE) {
  mm <- as.matrix(mat)
  out <- t(mm) %*% mm

  if(!weighted) {
    diag(out) <- 0
  }

  if(!weighted) {
    out[out > 0] <- 1
  }

  out
}
