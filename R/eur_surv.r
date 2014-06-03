library(survival)
library(plyr)
library(igraph)

source('../R/paleo_surv.r')

source('../R/cosmo_prov.r')
source('../R/oxygen_curve.r')

erdur <- read.csv('../data/euro-ranges.csv', stringsAsFactors = FALSE)
names(erdur) <- c('genus', 'species', 'fad', 'lad', 
                  'collections', 'abundance', 'geo.mean.ab')

load('../data/body_mass.rdata')  # body mass data

bi <- with(erdur, binom.make(genus, species))
erdur <- erdur[bi %in% eur$name.bi, ]
erdur$name.bi <- with(erdur, binom.make(genus, species))

ecol <- cbind(data.frame(taxa = eur$name.bi, stringsAsFactors = FALSE),
              diet = eur$comdiet, move = eur$comlife)

ecol <- ecol[order(ecol$taxa), ]
er.ecol <- ecol[!duplicated(ecol$taxa), ]

#er.ecol <- er.ecol[er.ecol$taxa %in% er.mass$name, ]
#er.ecol <- cbind(er.ecol, mass = er.mass$value)
#erdur <- erdur[erdur$name.bi %in% er.mass$name, ]
#gg <- er.mass$value != 0
#er.ecol <- er.ecol[gg, ]
#erdur <- erdur[gg, ]

# occupancy
taxa.occ <- lapply(eurwin, function(x) {
                   occupancy(x, membership = membership(infomap.community(x)))})
taxa.occ <- Reduce(rbind, taxa.occ)
rewrite <- order(as.character(taxa.occ$taxa))
taxa.occ <- taxa.occ[rewrite, ]
sp.occ <- split(taxa.occ, taxa.occ$taxa)
mean.occ <- melt(lapply(sp.occ, function(x) mean(x[, 1])))
names(mean.occ) <- c('occ', 'taxa')
mean.occ <- mean.occ[order(mean.occ$taxa), ]
cv.occ <- melt(lapply(sp.occ, function(x) var(x[, 1]) / mean(x[, 1])))
names(cv.occ) <- c('cv.occ', 'taxa')
cv.occ <- cv.occ[order(cv.occ$taxa), ]

# climate
isotope.match <- Map(getclimate, erdur$fad, erdur$lad)
names(isotope.match) <- erdur$name.bi
mean.climate <- melt(lapply(isotope.match, function(x) mean(x, na.rm = TRUE)))
names(mean.climate) <- c('climate', 'taxa')
cv.climate <- melt(lapply(isotope.match, function(x) {
                           var(x, na.rm = TRUE) / abs(mean(x, na.rm = TRUE))}))
names(cv.climate) <- c('cv.climate', 'taxa')

# exclude taxa that originate after cutoff
young <- which(nadur[, 3] <= 2)
nadur <- nadur[-young, ]
na.ecol <- na.ecol[-young, ]

# remove 0 survival times
rms <- which(abs(nadur[, 3] - nadur[, 4]) <= 0.1)
nadur <- nadur[-rms, ]
na.ecol <- na.ecol[-rms, ]

# exclude taxa that origierte after cutoff
young <- which(erdur[, 3] <= 2)
erdur <- erdur[-young, ]
er.ecol <- er.ecol[-young, ]

# remove 0 survival times
rms <- which(abs(erdur[, 3] - erdur[, 4]) <= 0.1)
erdur <- erdur[-rms, ]
er.ecol <- er.ecol[-rms, ]

er.surv <- paleosurv(fad = erdur[, 3], lad = erdur[, 4],
                     start = 66, end = 2)


# change to generic level
ergen <- ddply(erdur, .(genus), summarize,
               fad = max(fad),
               lad = min(lad))
ergen.surv <- paleosurv(fad = ergen[, 2], lad = ergen[, 3], start = 66, end = 2)


# occupancy
taxa.occ <- lapply(ergenwin, function(x) {
                   occupancy(x, membership = membership(infomap.community(x)))})
taxa.occ <- Reduce(rbind, taxa.occ)
rewrite <- order(as.character(taxa.occ$taxa))
taxa.occ <- taxa.occ[rewrite, ]
sp.occ <- split(taxa.occ, taxa.occ$taxa)
mean.occ <- melt(lapply(sp.occ, function(x) mean(x[, 1])))
names(mean.occ) <- c('occ', 'taxa')
mean.occ <- mean.occ[order(mean.occ$taxa), ]
cv.occ <- melt(lapply(sp.occ, function(x) var(x[, 1]) / mean(x[, 1])))
names(cv.occ) <- c('cv.occ', 'taxa')

# climate
generic.isotope <- Map(getclimate, ergen$fad, ergen$lad)
names(generic.isotope) <- ergen$genus
gmean.climate <- melt(lapply(generic.isotope, function(x) mean(x, na.rm = TRUE)))
names(gmean.climate) <- c('climate', 'taxa')
gcv.climate <- melt(lapply(generic.isotope, function(x) {
                            var(x, na.rm = TRUE) / abs(mean(x, na.rm = TRUE))}))
names(gcv.climate) <- c('cv.climate', 'taxa')


er.genecol <- cbind(er.ecol, genus = erdur$genus)
er.genecol <- ddply(er.genecol, .(genus), summarize,
                    diet = names(which.max(table(diet))),
                    move = names(which.max(table(move))))

