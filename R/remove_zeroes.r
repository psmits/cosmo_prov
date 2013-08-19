#' Remove sites or taxa with zero samples
#'
#' This functions takes a multivariate abundance matrix and removes any and all
#' sites or taxa that have samples of zero. This tolerance level is user 
#' adjustable so that, for example, singleton taxa could be removed.
#'
#' TODO find a matrix that actually *needs* this function
#'
#' @param pbdb object of class pbdb
#' @param tol tolerance level for removal
#' @keywords
#' @export
#' @examples
remove.low <- function(pbdb, tol = 1) {
  mat <- pbdb$occurence
  if(is.element('ecology', names(pbdb))) {
    eco <- pbdb$ecology
  } 
  if(is.element('geology', names(pbdb))) {
    geo <- pbdb$geology
  }

  while(any(colSums(mat) <= tol) || any(rowSums(mat) <= tol)) {
    cc <- colSums(mat) <= tol
    rr <- rowSums(mat) <= tol

    mat <- mat[!rr, !cc]
    eco <- eco[!rr, ]
    geo <- geo[, !cc]
  }

  pbdb$occurence <- mat
  pbdb$ecology <- eco
  pbdb$geology <- geo

  pbdb
}


