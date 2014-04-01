library(survival)
library(plyr)

source('../R/paleo_surv.r')

source('../R/na_mung.r')

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

best <- na.ecol[na.ecol$taxa %in% na.mass$name, ]
best <- cbind(best, mass = na.mass$value)
bestdur <- nadur[nadur$name.bi %in% na.mass$name, ]

# exclude taxa that originate after cutoff
young <- which(nadur[, 3] <= 2)
nadur <- nadur[-young, ]
na.ecol <- na.ecol[-young, ]
byong <- which(bestdur[, 3] <= 2)
bestdur <- bestdur[-byong, ]
best <- best[-byong, ]

# remove 0 survival times
rms <- which(abs(nadur[, 3] - nadur[, 4]) <= 0.1)
nadur <- nadur[-rms, ]
na.ecol <- na.ecol[-rms, ]
rms <- which(abs(bestdur[, 3] - bestdur[, 4]) <= 0.1)
bestdur <- bestdur[-rms, ]
best <- best[-rms, ]

na.surv <- paleosurv(fad = nadur[, 3], lad = nadur[, 4],
                     start = 66, end = 2)
best.surv <- paleosurv(fad = bestdur[, 3], lad = bestdur[, 4],
                       start = 66, end = 2)


# change to generic level
nagen <- ddply(nadur, .(genus), summarize,
               fad = max(fad),
               lad = min(lad))
nagen.surv <- paleosurv(fad = nagen[, 2], lad = nagen[, 3], start = 66, end = 2)

na.genecol <- ddply(dat, .(occurrence.genus_name), summarize,
                    diet = names(which.max(table(comdiet))),
                    move = names(which.max(table(comlife))))

na.genecol <- na.genecol[na.genecol[, 1] %in% nagen$genus, ]
