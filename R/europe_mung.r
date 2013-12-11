library(plyr)
library(reshape2)
library(taxize)

source('../R/clean_pbdb.r')
source('../R/mung_help.r')
source('../R/collapse_names.r')
source('../R/remove_zeroes.r')
source('../R/taxon_names.r')
source('../R/my_get_eolid.r')
eol.key = '2a9932f264f3f0421db36158b6e785b535c6da0e'

eur <- read.csv('../data/euro-occs.csv', stringsAsFactors = FALSE)

# remove specimens that don't have time assigned
eur <- eur[!is.na(eur$ma_mid), ]
eur <- eur[!is.na(eur$ma_max), ]
eur <- eur[!is.na(eur$ma_min), ]

# make sure i'm entirely in the Cenozoic via max
kpg <- 65.6 
eur <- eur[!eur$ma_max > kpg, ]

# lat/long fix
eur$paleolatdec <- as.numeric(eur$paleolatdec)
eur$paleolngdec <- as.numeric(eur$paleolngdec)
eur <- eur[!(is.na(eur$paleolatdec) | is.na(eur$paleolngdec)), ]

# remove all the sp.-s
#grep('sp', x = dat$occurence.species_name, perl = TRUE)
eur <- eur[eur$occurrence.species_name != 'sp.', ]  # change to a grep 

binm <- with(eur, binom.make(occurrence.genus_name, occurrence.species_name))
eur <- cbind(eur, name.bi = binm)
eur$name.bi <- as.character(eur$name.bi)

# exclude aquatic and volant nonsense
aq <- c('Cetacea', 'Desmostylia', 'Sirenia', 'Chiroptera')
eur <- eur[!(eur$order_name %in% aq), ]
lf <- c('amphibious', 'volant', 'aquatic', 'glinding')
eur <- eur[!(eur$life_habit %in% lf), ]

# diet assignments
herb <- c('herbivore', 'grazer', 'browser', 'folivore', 'granivore')
omm <- c('frugivore', 'omnivore')
car <- c('carnivore', 'insectivore')
eur$comdiet <- eur$diet1
eur$comdiet[eur$diet1 %in% herb] <- 'herb'
eur$comdiet[eur$diet1 %in% omm] <- 'omni'
eur$comdiet[eur$diet1 %in% car] <- 'carni'

# locomotor assignments
tree <- c('arboreal')
ground <- c('ground dwelling', 'semifossorial', 'fossorial', 'saltatorial')
eur$comlife <- eur$life_habit
eur$comlife[eur$life_habit %in% tree] <- 'arboreal'
eur$comlife[eur$life_habit %in% ground] <- 'ground dwelling'

# assign every occurence to a 2 My bin
bins <- seq(from = 0, to = 66, by = 2)
bins <- cbind(top = bins[-1], bot = bins[-length(bins)])
eur$bins <- rep(NA, nrow(eur))
for (ii in seq(nrow(bins))) {
  out <- which(eur$ma_mid < bins[ii, 1] & eur$ma_mid >= bins[ii, 2])
  eur$bins[out] <- bins[ii, 1]
}

# 2x2, 5x5, 10x10 
eur$gid <- with(eur, grid.id(paleolatdec, paleolngdec, 2))

# remove duplicates at grid locations in each bin
db <- split(eur, eur$bins)
dbg <- lapply(db, function(x) split(x, x$gid))
uu <- lapply(dbg, function(x) {
             lapply(x, function(y) {
                    dup <- duplicated(y$name.bi)
                    y[!dup, ]})})
uu <- lapply(uu, function(x) {
             rms <- lapply(x, nrow) == 0
             x[!rms]})
uu <- lapply(uu, function(x) Reduce(rbind, x))
eur <- Reduce(rbind, uu)

# taxonomic information
# get for all unique genera
gen <- unique(eur$occurrence.genus_name)
un <- unique(eur$name.bi)

gen <- gen[-which(gen %in% c('Vulpes', 'Myxomygale'))]

new.tax <- grab.heir(gen, key = eol.key)
new.tax <- apply(new.tax, 2, function(x) {
                 gsub(pattern = '\\s(.*)', x = x, 
                      perl = TRUE, replacement = '')})
bad.orders <- unique(new.tax[, 1])[c(10, 14:length(unique(new.tax[, 1])))]
new.tax <- new.tax[!(new.tax[, 1] %in% bad.orders), ]

for(ii in seq(nrow(new.tax))) {
  rp <- which(eur$occurrence.genus_name == new.tax[ii, 3])
  eur[rp, 'order_name'] <- new.tax[ii, 1]
  eur[rp, 'family_name'] <- new.tax[ii, 2]
}
