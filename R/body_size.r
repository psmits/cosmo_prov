library(stringr)
library(plyr)

source('../R/taxon_names.r')
source('../R/na_mung.r')
source('../R/europe_mung.r')

pbdb <- read.csv('../data/mam_avg_measure.txt', stringsAsFactors = FALSE)

now <- read.csv('../data/now_database.csv', stringsAsFactors = FALSE)
pan <- read.delim('../data/pantheria.txt', stringsAsFactors = FALSE)
tomiya <- read.csv('../data/tomiya_mass.csv', stringsAsFactors = FALSE)

pbdb <- pbdb[pbdb$measurement == 'mass', ]

now$BODYMASS[!(grepl('[0-9]', now$BODYMASS))] <- NA
now$BODYMASS <- as.numeric(now$BODYMASS)

pan$mass <- pan$X5.5_AdultBodyMass_g_EXT
pan$mass[pan$mass < 0] <- NA

get.val <- function(data, target, val) {
  mat <- which(data %in% target)
  out <- data.frame(name = data[mat], value = val[mat])
  out <- unique(out[!is.na(out[, 2]), ])
  out
}

na.pbdb <- get.val(pbdb$species, dat$name.bi, as.numeric(pbdb$mean))
er.pbdb <- get.val(pbdb$species, eur$name.bi, as.numeric(pbdb$mean))

# NOW
nowbi <- binom.make(now$GENUS, now$SPECIES)
na.now <- get.val(nowbi, dat$name.bi, now$BODYMASS)
er.now <- get.val(nowbi, eur$name.bi, now$BODYMASS)

# Pantheria
panbi <- pan$MSW05_Binomial
na.pan <- get.val(panbi, dat$name.bi, pan$mass)
er.pan <- get.val(panbi, eur$name.bi, pan$mass)

# Susumu's work
tomiya$Taxon <- str_replace(tomiya$Taxon, '_', ' ')
na.tom <- get.val(tomiya$Taxon, dat$name.bi, exp(tomiya$LnMass))
er.tom <- get.val(tomiya$Taxon, eur$name.bi, exp(tomiya$LnMass))

na.mass <- rbind(na.pbdb, na.now, na.pan, na.tom)
er.mass <- rbind(na.pbdb, er.now, er.pan, er.tom)

# clear out duplicates
na.mass <- na.mass[!duplicated(na.mass[, 1]), ]
er.mass <- er.mass[!duplicated(er.mass[, 1]), ]

save(na.mass, er.mass,
     file = '../data/body_mass.rdata')
