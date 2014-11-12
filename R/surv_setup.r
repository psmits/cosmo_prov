library(survival)
library(plyr)
library(igraph)

source('../R/paleo_surv.r')

source('../R/cosmo_prov.r')
source('../R/oxygen_curve.r')

nadur <- read.csv('../data/mam-ranges.csv', stringsAsFactors = FALSE)
names(nadur) <- c('genus', 'species', 'fad', 'lad', 
                  'collections', 'abundance', 'geo.mean.ab')

load('../data/body_mass.rdata')  # body mass data

bi <- with(nadur, binom.make(genus, species))
nadur <- nadur[bi %in% dat$name.bi, ]
nadur$name.bi <- with(nadur, binom.make(genus, species))

ecol <- cbind(data.frame(taxa = dat$name.bi, stringsAsFactors = FALSE),
              diet = dat$comdiet, move = dat$comlife)
ecol <- ecol[order(ecol$taxa), ]
na.ecol <- ecol[!duplicated(ecol$taxa), ]

na.ecol <- na.ecol[na.ecol$taxa %in% na.mass$name, ]
na.ecol <- cbind(na.ecol, mass = na.mass$value)
nadur <- nadur[nadur$name.bi %in% na.mass$name, ]

# occupancy
taxa.occ <- lapply(taxawin, function(x) {
                   occupancy(x, membership = membership(infomap.community(x)))})
taxa.occ <- Reduce(rbind, taxa.occ)
rewrite <- order(as.character(taxa.occ$taxa))
na.taxa.occ <- taxa.occ <- taxa.occ[rewrite, ]
sp.occ <- split(taxa.occ, taxa.occ$taxa)
mean.occ <- melt(lapply(sp.occ, function(x) mean(x[, 1])))
names(mean.occ) <- c('occ', 'taxa')
na.mean.occ <- mean.occ[order(mean.occ$taxa), ]

cv.occ <- melt(lapply(sp.occ, function(x) var(x[, 1]) / mean(x[, 1])))
names(cv.occ) <- c('cv.occ', 'taxa')
na.cv.occ <- cv.occ[order(cv.occ$taxa), ]
# climate
isotope.match <- Map(getclimate, nadur$fad, nadur$lad)
names(isotope.match) <- nadur$name.bi
na.mean.climate <- melt(lapply(isotope.match, function(x) mean(x, na.rm = TRUE)))
names(na.mean.climate) <- c('climate', 'taxa')
na.cv.climate <- melt(lapply(isotope.match, function(x) {
                           var(x, na.rm = TRUE) / abs(mean(x, na.rm = TRUE))}))
names(na.cv.climate) <- c('cv.climate', 'taxa')
no.climate <- is.na(na.mean.climate[, 1])
nadur <- nadur[!no.climate, ]
na.ecol <- na.ecol[!no.climate, ]
na.ecol$climate <- na.mean.climate[!no.climate, 1]
na.ecol$cv.li <- na.cv.climate[!no.climate, 1]

# exclude taxa that originate after cutoff
bin.range <- ddply(dat, .(name.bi), summarize, 
                   old <- max(bins),
                   young <- min(bins))
bin.range <- bin.range[bin.range[, 1] %in% na.ecol$taxa,]

dur <- (bin.range[, 2] - bin.range[, 3]) / 2 + 1
ext <- as.numeric(bin.range[, 3] != 2)
occ <- na.mean.occ[na.mean.occ[, 2] %in% bin.range[, 1], 1]

cohort <- bin.range[, 2] / 2  # find the cohorts
