library(stringr)
library(plyr)
source('../R/predict_mass.r')

source('../R/taxon_names.r')
source('../R/na_mung.r')
source('../R/europe_mung.r')

pbdb <- read.csv('../data/mam_avg_measure.txt', stringsAsFactors = FALSE)

now <- read.csv('../data/now_database.csv', stringsAsFactors = FALSE)
pan <- read.delim('../data/pantheria.txt', stringsAsFactors = FALSE)
tomiya <- read.csv('../data/tomiya_mass.csv', stringsAsFactors = FALSE)

# play with the pbdb measures
#first.match <- str_detect(pbdb$position, 'm1')
#problem.dash <- str_detect(pbdb$position, '\\-')
#good <- first.match & !problem.dash
#pbdb$position[first.match & !problem.dash]
#simple solution
tooth.spec <- pbdb$species[(pbdb$position == 'm1' | pbdb$position == 'Lower m1') 
                           & pbdb$measurement == 'area']
mass.spec <- pbdb$species[pbdb$measurement == 'mass']
update.mass <- tooth.spec[!(tooth.spec %in% mass.spec)]
update.block <- pbdb[pbdb$species %in% update.mass, ]
update.block <- update.block[update.block$measurement == 'area', ]
tooth.mean <- ddply(update.block, .(species), summarize,
                    mean(as.numeric(mean)))
tooth.mean$mass <- predmass(tooth.mean[, 2])

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

# my predictions from above
na.me <- get.val(tooth.mean[, 1], dat$name.bi, tooth.mean$mass)
er.me <- get.val(tooth.mean[, 1], eur$name.bi, tooth.mean$mass)

na.mass <- rbind(na.pbdb, na.now, na.pan, na.tom, na.me)
er.mass <- rbind(er.pbdb, er.now, er.pan, er.tom, er.me)


# clear out duplicates
na.mass <- na.mass[!duplicated(na.mass[, 1]), ]
er.mass <- er.mass[!duplicated(er.mass[, 1]), ]

save(na.mass, er.mass,
     file = '../data/body_mass.rdata')
