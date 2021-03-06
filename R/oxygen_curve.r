zachos <- read.csv('../data/2008_zachos_data.csv', header = T)
tits <- c('site', 'age', 'genus', 'd18o', 'd13c', 'd18o.5', 'd13c.5')
zachos <- zachos[, 1:7]
names(zachos) <- tits

# time bin means
# is there an easier way to do this with ddply?
bin <- seq(from = 0, to = 66, by = 2)
oo <- array(dim = length(bin) - 1)
for(ii in seq(length(bin) - 1)) {
  ww <- zachos$age > bin[ii] & zachos$age <= bin[ii + 1]
  oo[ii] <- mean(zachos$d18o[ww], na.rm = TRUE)
}

names(oo) <- (bin + 2)[seq(length(bin) - 1)]

#' Get climate proxy for duration of a taxon
#'
#' Given the Zachos et al 2008 dO^18 get all isotope values from a taxon's duration.
#'
#' @param fad numeric first apperance date. Now is 0.
#' @param lad numeric last apperance date. Now is 0.
#' @export
#' @return vector of isotope values
getclimate <- function(fad, lad) {
  times <- zachos$age

  isotope <- zachos$d18o[times <= fad & times >= lad]
  isotope
}
