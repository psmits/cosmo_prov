library(stringr)
library(xtable)
library(plyr)
source('../R/predict_mass.r')

source('../R/taxon_names.r')
source('../R/na_mung.r')
source('../R/europe_mung.r')

load('../data/update_taxonomy.rdata')

pbdb <- read.csv('../data/mam_avg_measure.txt', stringsAsFactors = FALSE)
now <- read.csv('../data/now_database.csv', stringsAsFactors = FALSE)
pan <- read.delim('../data/pantheria.txt', stringsAsFactors = FALSE)
tomiya <- read.csv('../data/tomiya_mass.csv', stringsAsFactors = FALSE)
smith <- read.csv('../data/smith2003.csv', stringsAsFactors = FALSE)
brook <- read.csv('../data/brook_bowman_2004.csv', stringsAsFactors = FALSE)
raia <- read.csv('../data/raia_geb_mass.csv', stringsAsFactors = FALSE)

smits <- read.csv('../data/found_mass.csv', stringsAsFactors = FALSE)

get.val <- function(data, target, val) {
  mat <- which(data %in% target)
  out <- data.frame(name = data[mat], value = val[mat])
  out <- unique(out[!is.na(out[, 2]), ])
  out
}

# what do i have estimates for already?
# PBDB
gpbdb <- pbdb[pbdb$measurement == 'mass', ]
na.pbdb <- get.val(gpbdb$species, dat$name.bi, as.numeric(gpbdb$mean))
er.pbdb <- get.val(gpbdb$species, eur$name.bi, as.numeric(gpbdb$mean))
# NOW
now$BODYMASS[!(grepl('[0-9]', now$BODYMASS))] <- NA
now$BODYMASS <- as.numeric(now$BODYMASS)
nowbi <- binom.make(now$GENUS, now$SPECIES)
na.now <- get.val(nowbi, dat$name.bi, now$BODYMASS)
er.now <- get.val(nowbi, eur$name.bi, now$BODYMASS)
# Pantheria
pan$mass <- pan$X5.5_AdultBodyMass_g_EXT
pan$mass[pan$mass < 0] <- NA
panbi <- pan$MSW05_Binomial
na.pan <- get.val(panbi, dat$name.bi, pan$mass)
er.pan <- get.val(panbi, eur$name.bi, pan$mass)
# Smith et al. 2003
smithbi <- binom.make(smith$Genus, smith$Species)
na.smith <- get.val(smithbi, dat$name.bi, exp(smith$LogMass))
er.smith <- get.val(smithbi, eur$name.bi, exp(smith$LogMass))
# brook and bownman 2004
brook$mass <- 10^brook[, 3]
na.brook <- get.val(brook$Name, dat$name.bi, brook$mass)
er.brook <- get.val(brook$Name, eur$name.bi, brook$mass)
# Susumu's work
tomiya$Taxon <- str_replace(tomiya$Taxon, '_', ' ')
na.tom <- get.val(tomiya$Taxon, dat$name.bi, exp(tomiya$LnMass))
er.tom <- get.val(tomiya$Taxon, eur$name.bi, exp(tomiya$LnMass))
# raia et al. 2010
raia$mass <- 10^raia[, 3]
na.raia <- get.val(raia$Species, dat$name.bi, raia$mass)
er.raia <- get.val(raia$Species, eur$name.bi, raia$mass)
# mass values i've found
smits$value <- as.numeric(smits$value)
me.good <- smits[smits$measure == 'mass', ]
na.me <- get.val(me.good$species, dat$name.bi, me.good$value)
er.me <- get.val(me.good$species, eur$name.bi, me.good$value)

# clear out duplicates
na.mass <- rbind(na.pbdb, na.now, na.smith, na.brook, na.pan, na.tom, na.raia, na.me)
er.mass <- rbind(er.pbdb, er.now, er.smith, er.brook, er.pan, er.tom, er.raia, er.me)
na.mass <- na.mass[!duplicated(na.mass[, 1]), ]
er.mass <- er.mass[!duplicated(er.mass[, 1]), ]

north.source <- list(na.pbdb, na.now, na.smith, na.brook, na.pan, na.tom, 
                     na.raia, na.me)
sources <- list('PBDB', 'NOW', 'Smith2004', 'Brook2004', 'PanTheria', 
                'Tomiya2013', 'Raia2010', 'this study')
north.source <- Map(function(x, y) cbind(x, source = rep(y, nrow(x))), 
                    north.source, sources)
north.source <- Reduce(rbind, north.source)
north.source <- north.source[!(duplicated(north.source[, 1])), ]


# for those remaining
missing.mass <- unique(dat$name.bi[!(dat$name.bi %in% na.mass$name)],
                       eur$name.bi[!(eur$name.bi %in% er.mass$name)])
# remove the bad names
me.bad <- smits$species[which(smits$notes != '' & 
                              smits$notes != 'alt')]
missing.mass <- missing.mass[-which(missing.mass %in% me.bad)]
# update the taxonomy
na.uptax <- na.tax[na.tax$name.bi %in% missing.mass, ]
er.uptax <- er.tax[er.tax$name.bi %in% missing.mass, ]
# split the data by order
ungulates <- c('Proboscidea', 'Perissodactyla', 'Cetartiodactyla', 
               'Notounguluata', 'Artiodactyla', 'Condylartha', 'Dinocerata')
insectivore <- c('Eulipotyphla', 'Soricomorpha', 'Lipotyphla', 'Leptictida',
                 'Cimolesta', 'Afrosoricida', 'Trituberculata', 'Macroscelidea')
carnivore <- c('Carnivora', 'Credonta', 'Mesonychia', 'Miacoidae')
lagomorph <- 'Lagomorpha'
rodentia <- 'Rodentia'
marsupial <- 'Didelphimorphia'

ords <- c(na.uptax$order_name, er.uptax$order_name)
groups <- rep('general', length(ords))
groups[ords %in% ungulates] <- 'ungulates'
groups[ords %in% insectivore] <- 'insectivore'
groups[ords %in% carnivore] <- 'carnivore'
groups[ords %in% lagomorph] <- 'lagomorph'
groups[ords %in% rodentia] <- 'rodentia'
groups[ords %in% marsupial] <- 'marsupial'
ords <- split(c(na.uptax$name.bi, er.uptax$name.bi), groups)

# condense measurement data
# pbdb
good.position <- c('m1', 'Lower m1', 'Lower m1 ', 'M1', 'Upper M1', 'Upper M1 ',
                   'm2', 'lower M2', 'Lower m2', 'mandible', 'length of mandible',
                   'Mandibular', 'skull', 
                   'Entire skull length measured from I1 to occipital condyles')

good.measure <- pbdb[pbdb$position %in% good.position, ]
good.mean <- ddply(good.measure, .(species, position, measurement), summarize,
                   value = mean(as.numeric(mean)))
mandib <- which(good.mean$position %in% c('mandible', 'length of mandible',
                                          'Mandibular'))
rms <- mandib[which(good.mean$measurement[mandib] != 'length')]
if(length(rms) > 0) good.mean <- good.mean[-rms, ]
skull <- which(good.mean$position %in% c('skull', 'Entire skull length measured 
                                         from I1 to occipital condyles'))
rms <- skull[which(good.mean$measurement[skull] != 'length')]
if(length(rms) > 0) good.mean <- good.mean[-rms, ]
# mine
me.measure <- smits[smits$measure != '' & smits$measure != 'mass', ]
# combined
names(good.mean)[2:3] <- c('part', 'measure')
measures <- rbind(me.measure[, c(1, 3:5)], good.mean[, 1:4])

# predict mass
# ungulate
ung.mass <- ungulate.mass(ords$ungulates, measures)
# carnivores
car.mass <- carnivore.mass(ords$carnivore, measures)
# lagomorphs
lag.mass <- lagomorph.mass(ords$lagormoph, measures)
# insectivores
ins.mass <- insectivore.mass(ords$insectivore, measures)
# rodents
rod.mass <- rodentia.mass(ords$rodentia, measures)
# marsupials
mar.mass <- marsupial.mass(ords$marsupial, measures)
# general 
gen.mass <- general.mass(c(na.uptax$name.bi, er.uptax$name.bi),
                         measures)
est.mass <- rbind(ung.mass, car.mass, lag.mass, ins.mass, rod.mass, mar.mass, gen.mass)
est.mass <- est.mass[!duplicated(est.mass$species), ]
na.est <- get.val(est.mass$species, dat$name.bi, est.mass$value)
er.est <- get.val(est.mass$species, eur$name.bi, est.mass$value)
na.mass <- rbind(na.mass, na.est)
er.mass <- rbind(er.mass, er.est)

na.missing <- unique(dat$name.bi[!(dat$name.bi %in% na.mass$name)])
eur.missing <- unique(eur$name.bi[!(eur$name.bi %in% er.mass$name)])
missing.taxa <- sort(unique(c(eur.missing, na.missing)))
#write.csv(missing.taxa, '../data/unknown_mass.csv', row.names = FALSE)

save(na.mass, er.mass, file = '../data/body_mass.rdata')


# make the sources table
north.source <- rbind(north.source, 
                      cbind(na.est, source = rep('PBDB + regression', 
                                                 nrow(na.est))))
founds <- me.good$species %in% north.source[north.source$source == 
                                            'this study', 1]
new.source <- me.good[founds, c('species', 'source')]
north.source$source <- as.character(north.source$source)
north.source[north.source[, 1] %in% new.source[, 1], 
             'source'] <- new.source[, 2]
#unique(north.source$source)
the.fixer <- matrix(c('Smith', 'cite{Smith2004}', 
                      'Brook', 'cite{Brook2004a}', 
                      'Tomiya', 'cite{Tomiya2013}', 
                      'Raia', 'cite{Raia2012f}', 
                      'McKenna', 'cite{McKenna2011}', 
                      'Wilson et', 'cite{Wilson2012}', 
                      'Sorkin', 'cite{Sorkin2008}', 
                      'Conroy', 'cite{Controy1987}', 
                      'ADW', 'Animal Diversity Web', 
                      'Soligo', 'cite{Soligo2006}', 
                      'MacFadden', 'cite{MacFadden1986}', 
                      'McDonald', 'cite{McDonald1995}', 
                      'Martin', 'cite{Martin2002a}', 
                      'EOL', 'Encyclopedia of Life', 
                      'Strait', 'cite{Strait2001}', 
                      'Egi', 'cite{Egi2001}', 
                      'Torregrosa', 'cite{Torregrosa2010}'), 
                    ncol = 2, byrow = TRUE)
ss <- '(.*)$'
the.fixer[, 1] <- paste0(the.fixer[, 1], ss)

for(ii in seq(nrow(the.fixer))) {
  north.source$source <- str_replace(north.source$source, 
                                     the.fixer[ii, 1], the.fixer[ii, 2])
}
names(north.source) <- c('Species', 'Mass (g)', 'Source')
mass.table <- xtable(north.source, label = 'tab:mass_data')
print.xtable(mass.table, 
             file = '../doc/na_surv/mass_data.tex',
             include.rownames = FALSE, 
             sanitize.text.function = identity)
