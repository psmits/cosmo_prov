library(survival)
library(plyr)

source('../R/paleo_surv.r')

source('../R/europe_mung.r')

erdur <- read.csv('../data/euro-ranges.csv', stringsAsFactors = FALSE)
names(erdur) <- c('genus', 'species', 'fad', 'lad', 
                  'collections', 'abundance', 'geo.mean.ab')

bi <- with(erdur, binom.make(genus, species))
erdur <- erdur[bi %in% eur$name.bi, ]
erdur$name.bi <- with(erdur, binom.make(genus, species))

ecol <- cbind(data.frame(taxa = eur$name.bi, stringsAsFactors = FALSE),
              diet = eur$comdiet, move = eur$comlife)

ecol <- ecol[order(ecol$taxa), ]
er.ecol <- ecol[!duplicated(ecol$taxa), ]

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

genecol <- cbind(data.frame(genus = eur$occurrence.genus_name,
                            stringsAsFactors = FALSE),
                 diet = eur$comdiet, move = eur$comlife)
genecol <- genecol[order(genecol$genus), ]
er.genecol <- genecol[!duplicated(genecol$genus), ]

er.genecol <- er.genecol[er.genecol$genus %in% ergen$genus, ]
