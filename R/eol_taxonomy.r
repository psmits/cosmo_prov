#' Grab heirarchy
#'
#' @param taxon
#' @param key character; EoL api key
#' @return hierarchical information for taxa
#' @export
grab.heir <- function(taxon, key) {
  eol.id <- my.get_eolid(taxon, key = key)
  tax <- classification(eol.id, key = key)
  names(tax) <- taxon

  oandf <- lapply(tax, function(x) {
                  if(class(x) != 'character' & !is.na(x)){
                    if(any(x$rank == 'order', na.rm = TRUE)) {
                       pick <- which(x$rank == 'order')
                    } else {
                      pick <- which(x$rank == 'superorder') + 1
                    }
                    x[seq(pick, nrow(x)), ]}})
  oandf <- oandf[!laply(oandf, is.null)]

  oandf <- lapply(oandf, function(x) {
                  data.frame(lapply(x, as.character), 
                             stringsAsFactors = FALSE)})
  o <- list()
  for(ii in seq(length(oandf))) {
    o[[ii]] <- rbind(oandf[[ii]], c(names(oandf)[ii], 'genus'))
  }
  names(o) <- names(oandf)

  o <- lapply(o, function(x) {
              x <- t(x)
              colnames(x) <- x[2, ]
              x <- x[1, ]
              x})
  # exclude any not order, family or genus
  o <- lapply(o, function(x) {
              if(is.na(names(x)[1])) {
                names(x)[1] <- 'order'
                x
              } else { 
                x
              }})
  o <- lapply(o, function(x) {
              tt <- names(x) %in% c('order', 'family', 'genus')
              x[tt]})

  o <- Reduce(rbind, o)
  rownames(o) <- NULL
  o <- apply(o, 2, function(x) {
             gsub(pattern = '\\s(.*)', x = x,
                  perl = TRUE, replacement = '')})
  o
}
