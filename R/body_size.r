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
#na.smith <- get.val(smithbi, dat$name.bi, exp(smith$LogMass)) # old
#er.smith <- get.val(smithbi, eur$name.bi, exp(smith$LogMass)) # old
na.smith <- get.val(smithbi, dat$name.bi, 10^(smith$LogMass)) # correction
er.smith <- get.val(smithbi, eur$name.bi, 10^(smith$LogMass)) # correction
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

#oo <- measures[laply(str_split(measures[, 1], ' '), function(x) x[1]) %in% 'Acritoparamys', ]
#rodentia.mass(oo[, 1], oo)
# predict mass
# ungulate
ung.mass <- ungulate.mass(ords$ungulates, measures)
miohippus <- ords$ungulates[laply(str_split(ords$ungulates, ' '), 
                                  function(x) x[1]) == 'Miohippus']
ungulate.mass(miohippus, measures)

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
na.est <- get.val(est.mass$species, dat$name.bi, est.mass$mass) # old
er.est <- get.val(est.mass$species, eur$name.bi, est.mass$mass) # correction
#na.est <- get.val(est.mass$species, dat$name.bi, est.mass$value) # old
#er.est <- get.val(est.mass$species, eur$name.bi, est.mass$mass) # old
na.mass <- rbind(na.mass, na.est)
er.mass <- rbind(er.mass, er.est)

na.missing <- unique(dat$name.bi[!(dat$name.bi %in% na.mass$name)])
eur.missing <- unique(eur$name.bi[!(eur$name.bi %in% er.mass$name)])
missing.taxa <- sort(unique(c(eur.missing, na.missing)))
#write.csv(missing.taxa, '../data/unknown_mass.csv', row.names = FALSE)

save(na.mass, er.mass, file = '../data/body_mass.rdata')


# make the sources table
north.source$source <- as.character(north.source$source)

me.happy <- smits[smits$species %in% na.me$name, ]
north.source[north.source$name %in% me.happy$species, 'source'] <- 
  me.happy$source

north.source <- rbind(north.source, 
                      cbind(na.est, source = rep('PBDB + regression', 
                                                 nrow(na.est))))
founds <- me.measure[, 1] %in% north.source[, 1]
new.source <- me.measure[founds, c('species', 'source')]
north.source[north.source[, 1] %in% new.source[, 1], 
             'source'] <- new.source[, 2]

weird <- north.source[north.source$source == 'this study', ]
as.character(weird[, 1]) %in% me.measure$species


the.fixer <- matrix(c('Smith', 'cite{Smith2004}', 
                      'Brook', 'cite{Brook2004a}', 
                      'Tomiya', 'cite{Tomiya2013}', 
                      'Raia', 'cite{Raia2012f}', 
                      'McKenna', 'cite{McKenna2011}', 
                      'Wilson et', 'cite{Wilson2012}', 
                      'Sorkin', 'cite{Sorkin2008}', 
                      'Conroy', 'cite{Conroy1987}', 
                      'ADW', 'Animal Diversity Web', 
                      'Soligo', 'cite{Soligo2006}', 
                      'MacFadden', 'cite{MacFadden1986}', 
                      'McDonald', 'cite{McDonald1995}', 
                      'Martin', 'cite{Martin2002a}', 
                      'EOL', 'Encyclopedia of Life', 
                      'Strait', 'cite{Strait2001}', 
                      'Egi', 'cite{Egi2001}', 
                      'Torregrosa', 'cite{Torregrosa2010}',
                      'Osborn 1933', 'cite{Osborn1933}',
                      'Bloch et al', 'cite{Bloch2007}',
                      'Ferrusquia', 'cite{Ferrusquia-Villafranca2006}',
                      'Scott and', 'cite{Scott1940}',
                      'Scott 2003 J', 'cite{Scott2003a}',
                      'Skinner', 'cite{Skinner1972}',
                      'Williamson and Brusatte 2013', 'cite{Williamson2013}',
                      'Lofgren and Anad', 'cite{Lofgren2011}',
                      'Czaplewski', 'cite{Czaplewski2012}',
                      'Beatty and', 'cite{Beatty2009}',
                      'Van Valkenburgh 2007', 'cite{VanValkenburgh2007a}',
                      'Loomis 1911', 'cite{Loomis1911}',
                      'Wang 1994', 'cite{Wang1994a}',
                      'Carraway', 'cite{Carraway2010}', 
                      'Emry et al', 'cite{Emry2005}',
                      'Coombs 1979', 'cite{Coombs1979}',
                      'Mora and Zamora', 'cite{Mora2005}',
                      'Dawson and Beard', 'cite{Dawson2007}',
                      'Becker and White', 'cite{Becker1981}',
                      'Mellett 1969', 'cite{Mellett1969}',
                      'Lim et al', 'cite{Lim2001}',
                      'Fox and Scott', 'cite{Fox2011b}',
                      'Zakrzewski 1991', 'cite{Zakrzewski1991a}',
                      'Rose and Krause', 'cite{Rose1982a}',
                      'Wood 1962', 'cite{Wood1962}',
                      'Cope 1871', 'cite{Cope1871}',
                      'Jepsen 1932', 'cite{Jepsen1932}',
                      'Taylor and Webb', 'cite{Taylor1976}',
                      'Korth 1993 T', 'cite{Korth1993}',
                      'Chester and Beard', 'cite{Chester2012}',
                      'Madden and Storer 1985', 'cite{Madden1985}',
                      'Grohe', 'cite{Grohe2010}',
                      'Ivy 1990 C', 'cite{Ivy1990}',
                      'Patton and Taylor', 'cite{Patton1973}',
                      'Bever 2003', 'cite{Bever2003}',
                      'Mac Intyre 1966', 'cite{MacIntyre1966}',
                      'Stock 1948', 'cite{Stock1948}',
                      'Baskin 2004', 'cite{Baskin2004}',
                      'Kelly and Wood', 'cite{Kelley1954}',
                      'Rose et al 2011', 'cite{Rose2011a}',
                      'Johansen 1996', 'cite{Johanson1996}',
                      'Hay 1916', 'cite{Hay1916}',
                      'Macdonald 1951', 'cite{Macdonald1951}',
                      'Secord 2008', 'cite{Secord2008a}',
                      'Rich 1981', 'cite{Rich1981}',
                      'Stock 1937', 'cite{Stock1937}',
                      'Mihlbacher and Demere', 'cite{Mihlbachler2010}',
                      'Carranza-Castaneda and Walt', 'cite{Carranza-Castaneda1992}',
                      'Scott et al 2013', 'cite{Scott2013}',
                      'Tseng et al 2009', 'cite{Tseng2009}',
                      'Zack et al', 'cite{Zack2005}',
                      'Loomis 1932', 'cite{Loomis1932}',
                      'Scott et al 1937', 'cite{Scott1937}',
                      'Tedford et al 1994', 'cite{Tedford1994}',
                      'Wang et al 2014', 'cite{Wang2014}',
                      'Bjork 1970', 'cite{Bjork1970}',
                      'Korth 2010 P', 'cite{Korth2010}',
                      'Silcox and Willamson', 'cite{Silcox2012}',
                      'Zonneveld and Gu', 'cite{Zonneveld2003}', 
                      'Novacek 1977', 'cite{Novacek1977}',
                      'Scott 2004', 'cite{Scott2004}',
                      'McGrew 1939', 'cite{McGrew1939}',
                      'Mihlbacher and Solounias 2006', 'cite{Mihlbachler2006}',
                      'Cooke 2011', 'cite{Cooke2011}',
                      'Hall 1930', 'cite{Hall1930}',
                      'Hay 1969 P', 'cite{Hay1969}',
                      'Macdonald 1956', 'cite{Macdonald1956}',
                      'Clemens 2011', 'cite{Clemens2011}',
                      'Cassiliano 2008', 'cite{Cassiliano2008}',
                      'Dalquest 1978', 'cite{Dalquest1978}',
                      'Clemens and Williamson', 'cite{Clemens2005}',
                      'Rose et al 2013', 'cite{Rose2013a}',
                      'Kirk and Williams 2011', 'cite{Kirk2011}',
                      'Gidley 1920', 'cite{Gidley1920}',
                      'Robinson 1966', 'cite{Robinson1966}',
                      'Silcox et al 2008', 'cite{Silcox2008}',
                      'Sinclair 1915', 'cite{Sinclair1915}',
                      'Lillegraven 1977 Bull', 'cite{Lillegraven1977}',
                      'Albright 2000 Biostrat', 'cite{Albright2000}',
                      'White 1988', 'cite{White1988}',
                      'Wang et al 1999', 'cite{Wang1999}',
                      'Matthew 1901 B', 'cite{Matthew1901}',
                      'Worthman and Earle 18', 'cite{Wortman1893}',
                      'Gazin 1930 Thesis Carn', 'cite{Gazin1930}',
                      'Brown 1980 Transactions', 'cite{Brown1980}',
                      'Stirton 1932', 'cite{Stirton1932}',
                      'Simons 1960', 'cite{Simons1960}',
                      'Heissig 2012', 'cite{Heissig2012a}',
                      'Williamson et al 2012', 'cite{Williamson2012}',
                      'Dawson 2012 Swiss', 'cite{Dawson2012}',
                      'Baskin 2011 Pal', 'cite{Baskin2011}'
                      ), 
                    ncol = 2, byrow = TRUE)
ss <- '(.*)$'
the.fixer[, 1] <- paste0(the.fixer[, 1], ss)

for(ii in seq(nrow(the.fixer))) {
  north.source$source <- str_replace(north.source$source, 
                                     the.fixer[ii, 1], the.fixer[ii, 2])
}
names(north.source) <- c('Species', 'Mass (g)', 'Source')
north.source <- north.source[order(as.character(north.source[, 1])), ]

save(north.source, file = '../data/na_mass_table.rdata')
