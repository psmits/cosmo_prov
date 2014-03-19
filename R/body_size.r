library(stringr)
library(plyr)

source('../R/taxon_names.r')
source('../R/na_mung.r')
source('../R/europe_mung.r')

now <- read.csv('../data/now_database.csv', stringsAsFactors = FALSE)
pan <- read.delim('../data/pantheria.txt', stringsAsFactors = FALSE)

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

nowbi <- binom.make(now$GENUS, now$SPECIES)
na.now <- get.val(nowbi, dat$name.bi, now$BODYMASS)
er.now <- get.val(nowbi, eur$name.bi, now$BODYMASS)

panbi <- pan$MSW05_Binomial
na.pan <- get.val(panbi, dat$name.bi, pan$mass)
er.pan <- get.val(panbi, eur$name.bi, pan$mass)

na.mass <- rbind(na.now, na.pan)
er.mass <- rbind(er.now, er.pan)

save(na.mass, er.mass,
     file = '../data/body_mass.rdata')
