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
smith <- read.csv('../data/smith2003.csv', stringsAsFactors = FALSE)

brook <- read.csv('../data/brook_bowman_2004.csv', stringsAsFactors = FALSE)
brook$mass <- 10^brook[, 3]

raia <- read.csv('../data/raia_geb_mass.csv', stringsAsFactors = FALSE)
raia$mass <- 10^raia[, 3]

# play with the pbdb measures
#simple solution
# m1
tooth.spec <- pbdb$species[(pbdb$position == 'm1' | 
                            pbdb$position == 'Lower m1' |
                            pbdb$position == 'Lower m1 ') 
                           & pbdb$measurement == 'area']
mass.spec <- pbdb$species[pbdb$measurement == 'mass']
update.mass <- tooth.spec[!(tooth.spec %in% mass.spec)]
update.block <- pbdb[pbdb$species %in% update.mass, ]
update.block <- update.block[update.block$measurement == 'area', ]
tooth.mean <- ddply(update.block, .(species), summarize,
                    mean(as.numeric(mean)))
tooth.mean$mass <- predmass(tooth.mean[, 2])

# M1
upperM1 <- pbdb[(pbdb$position == 'M1' | 
                 pbdb$position == 'Upper M1' |
                 pbdb$position == 'Upper M1 '), ]
M1area <- which(upperM1$measurement == 'area')
M1area.mass <- massM1(as.numeric(upperM1$mean[M1area]))
M1len <- which(upperM1$measurement == 'length')
M1len.mass <- massM1len(as.numeric(upperM1$mean[M1len]))
M1.mass <- cbind(species = c(upperM1$species[M1area], upperM1$species[M1len]), 
                 mass = c(M1area.mass, M1len.mass))
M1.mass <- M1.mass[!duplicated(M1.mass[, 1]), ]

# m2
lowerm2 <- pbdb[(pbdb$position == 'm2' | 
                 pbdb$position == 'lower M2' | 
                 pbdb$position == 'Lower m2'), ]
m2len <- which(lowerm2$measurement == 'length')
m2len.mass <- massm2len(as.numeric(lowerm2$mean[m2len]))
m2.mass <- cbind(species = lowerm2$species[m2len],
                 mass = c(m2len.mass))
m2.mass <- m2.mass[!duplicated(m2.mass[, 1]), ]

# mandible length
mandible <- pbdb[pbdb$position %in% c('mandible', 'length of mandible', 'Mandibular'), ]
mandlen <- which(mandible$measurement == 'length')
mllen.mass <- massML(as.numeric(mandible$mean[mandlen]))
ml.mass <- cbind(species = mandible$species[mandlen],
                 mass = c(mllen.mass))
ml.mass <- ml.mass[!duplicated(ml.mass[, 1]), ]

# skull length
skull <- pbdb[pbdb$position %in% c('skull', 
                                   'Entire skull length measured from 
                                   I1 to occipital condyles'), ]
sklen <- which(skull$measurement == 'length')
sklen.mass <- massSL(as.numeric(skull$mean[sklen]))
sk.mass <- cbind(species = skull$species[sklen],
                 mass = c(sklen.mass))
sk.mass <- sk.mass[!duplicated(sk.mass[, 1]), ]


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

# Smith et al. 2003
smithbi <- binom.make(smith$Genus, smith$Species)
na.smith <- get.val(smithbi, dat$name.bi, exp(smith$LogMass))
er.smith <- get.val(smithbi, eur$name.bi, exp(smith$LogMass))

# brook and bownman 2004
na.brook <- get.val(brook$Name, dat$name.bi, brook$mass)
er.brook <- get.val(brook$Name, eur$name.bi, brook$mass)

# Susumu's work
tomiya$Taxon <- str_replace(tomiya$Taxon, '_', ' ')
na.tom <- get.val(tomiya$Taxon, dat$name.bi, exp(tomiya$LnMass))
er.tom <- get.val(tomiya$Taxon, eur$name.bi, exp(tomiya$LnMass))

# raia et al. 2010
na.raia <- get.val(raia$Species, dat$name.bi, raia$mass)
er.raia <- get.val(raia$Species, eur$name.bi, raia$mass)

# my predictions from above
na.me <- get.val(tooth.mean[, 1], dat$name.bi, tooth.mean$mass)
er.me <- get.val(tooth.mean[, 1], eur$name.bi, tooth.mean$mass)
na.me <- rbind(na.me, get.val(M1.mass[, 1], dat$name.bi, M1.mass[, 2]))
er.me <- rbind(er.me, get.val(M1.mass[, 1], eur$name.bi, M1.mass[, 2]))
na.me <- rbind(na.me, get.val(m2.mass[, 1], dat$name.bi, m2.mass[, 2]))
er.me <- rbind(er.me, get.val(m2.mass[, 1], eur$name.bi, m2.mass[, 2]))
na.me <- rbind(na.me, get.val(ml.mass[, 1], dat$name.bi, ml.mass[, 2]))
er.me <- rbind(er.me, get.val(ml.mass[, 1], eur$name.bi, ml.mass[, 2]))
na.me <- rbind(na.me, get.val(sk.mass[, 1], dat$name.bi, sk.mass[, 2]))
er.me <- rbind(er.me, get.val(sk.mass[, 1], eur$name.bi, sk.mass[, 2]))

na.me <- na.me[!duplicated(na.me[, 1]), ]
er.me <- er.me[!duplicated(er.me[, 1]), ]


na.mass <- rbind(na.pbdb, na.now, na.smith, na.brook, na.pan, na.tom, na.raia, na.me)
er.mass <- rbind(er.pbdb, er.now, er.smith, er.brook, er.pan, er.tom, er.raia, er.me)


# clear out duplicates
na.mass <- na.mass[!duplicated(na.mass[, 1]), ]
er.mass <- er.mass[!duplicated(er.mass[, 1]), ]


# which are missing
#missing.mass <- unique(dat$name[!(dat$name.bi %in% na.mass$name)],
#                       eur$name.bi[!(eur$name.bi %in% er.mass$name)])
#write.csv(missing.mass, file = '../data/unknown_mass.csv')
#oops <- read.csv('../data/found_mass.csv', stringsAsFactors = FALSE)

save(na.mass, er.mass,
     file = '../data/body_mass.rdata')
