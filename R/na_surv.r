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
taxa.occ <- taxa.occ[rewrite, ]
sp.occ <- split(taxa.occ, taxa.occ$taxa)
mean.occ <- melt(lapply(sp.occ, function(x) mean(x[, 1])))
names(mean.occ) <- c('occ', 'taxa')
mean.occ <- mean.occ[order(mean.occ$taxa), ]
cv.occ <- melt(lapply(sp.occ, function(x) var(x[, 1]) / mean(x[, 1])))
names(cv.occ) <- c('cv.occ', 'taxa')
cv.occ <- cv.occ[order(cv.occ$taxa), ]

# climate
isotope.match <- Map(getclimate, nadur$fad, nadur$lad)
names(isotope.match) <- nadur$name.bi
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

na.surv <- paleosurv(fad = nadur[, 3], lad = nadur[, 4],
                     start = 66, end = 2)


# change to generic level
nagen <- ddply(nadur, .(genus), summarize,
               fad = max(fad),
               lad = min(lad))
nagen.surv <- paleosurv(fad = nagen[, 2], lad = nagen[, 3], start = 66, end = 2)


generic.isotope <- Map(getclimate, nagen$fad, nagen$lad)
names(generic.isotope) <- nagen$genus
gmean.climate <- melt(lapply(generic.isotope, function(x) mean(x, na.rm = TRUE)))
names(gmean.climate) <- c('climate', 'taxa')
gcv.climate <- melt(lapply(generic.isotope, function(x) {
                            var(x, na.rm = TRUE) / abs(mean(x, na.rm = TRUE))}))
names(gcv.climate) <- c('cv.climate', 'taxa')


na.genecol <- cbind(na.ecol, genus = nadur$genus)
na.genecol <- ddply(na.genecol, .(genus), summarize,
                    diet = names(which.max(table(diet))),
                    move = names(which.max(table(move))),
                    mass = mean(mass))

