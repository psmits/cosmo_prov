library(plyr)
library(reshape2)
library(mapproj)
library(stringr)
library(dismo)
library(raster)
library(sp)
library(XML)
library(maptools)
library(foreign)
library(rgdal)

source('../R/clean_pbdb.r')
source('../R/mung_help.r')
source('../R/collapse_names.r')
source('../R/remove_zeroes.r')
source('../R/taxon_names.r')

dat <- read.csv('../data/mam-occs.csv', stringsAsFactors = FALSE)

# remove specimens that don't have time assigned
dat <- dat[!is.na(dat$ma_mid), ]
dat <- dat[!is.na(dat$ma_max), ]
dat <- dat[!is.na(dat$ma_min), ]

# make sure i'm entirely in the Cenozoic via max
kpg <- 65.6
dat <- dat[!dat$ma_max > kpg, ]

# lat/long fix
dat$paleolatdec <- as.numeric(dat$paleolatdec)
dat$paleolngdec <- as.numeric(dat$paleolngdec)
dat <- dat[!(is.na(dat$paleolatdec) | is.na(dat$paleolngdec)), ]

# remove all the sp.-s
#grep('sp', x = dat$occurence.species_name, perl = TRUE)
dat <- dat[dat$occurrence.species_name != 'sp.', ]  # change to a grep 
dat <- dat[is.na(str_match(dat$occurrence.species_name, '/')), ]

binm <- with(dat, binom.make(occurrence.genus_name, occurrence.species_name))
dat <- cbind(dat, name.bi = binm)
dat$name.bi <- as.character(dat$name.bi)

# exclude aquatic and volant nonsense
aq <- c('Cetacea', 'Desmostylia', 'Sirenia', 'Chiroptera')
seal <- c('Odobenidae', 'Otariidae', 'Phocidae', 'Desmatophocidae')
badgen <- c('Enaliarctos', 'Pacificotaria', 'Pinnarctidion', 'Pteronarctos', 'Wallia')
dat <- dat[!(dat$order_name %in% aq), ]
dat <- dat[!(dat$family_name %in% seal), ]
dat <- dat[!(dat$occurrence.genus_name %in% badgen), ]
lf <- c('amphibious', 'volant', 'aquatic', 'gliding')
dat <- dat[!(dat$life_habit %in% lf), ]

# diet assignments
herb <- c('herbivore', 'grazer', 'browser', 'folivore', 'granivore')
omm <- c('frugivore', 'omnivore')
car <- c('carnivore')#, 'insectivore')
insect <- c('insectivore')
dat$comdiet <- dat$diet1
dat$comdiet[dat$diet1 %in% herb] <- 'herb'
dat$comdiet[dat$diet1 %in% omm] <- 'omni'
dat$comdiet[dat$diet1 %in% car] <- 'carni'
dat$comdiet[dat$diet1 %in% insect] <- 'insect'

# locomotor assignments
tree <- c('arboreal')
ground <- c('ground dwelling', 'semifossorial', 'fossorial', 'saltatorial')
dat$comlife <- dat$life_habit
dat$comlife[dat$life_habit %in% tree] <- 'arboreal'
dat$comlife[dat$life_habit %in% ground] <- 'ground dwelling'

# assign every occurence to a 2 My bin
bins <- seq(from = 0, to = 66, by = 2)
bins <- cbind(top = bins[-1], bot = bins[-length(bins)])
dat$bins <- rep(NA, nrow(dat))
for (ii in seq(nrow(bins))) {
  out <- which(dat$ma_mid < bins[ii, 1] & dat$ma_mid >= bins[ii, 2])
  dat$bins[out] <- bins[ii, 1]
}

# spatial localitions
#wgs1984.proj <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
eq <- CRS("+proj=cea +lat_0=0 +lon_0=0 +lat_ts=30 +a=6371228.0 +units=km")
globe.map <- readShapeSpatial('../data/ne_10m_coastline.shp')  # from natural earth
proj4string(globe.map) <- eq

spatialref <- SpatialPoints(coords = dat[, c('paleolngdec', 'paleolatdec')],
                            proj4string = eq)  # wgs1984.proj

r <- raster(globe.map, nrows = 70, ncols = 34)
sp.ras <- trim(rasterize(spatialref, r))
membership <- cellFromXY(sp.ras, xy = dat[, c('paleolngdec', 'paleolatdec')])
dat$gid <- membership


#plot(sp.ras)
#plot(globe.map, add = TRUE)

# remove duplicates at grid locations in each bin
db <- split(dat, dat$bins)
dbg <- lapply(db, function(x) split(x, x$gid))
uu <- lapply(dbg, function(x) {
             lapply(x, function(y) {
                    dup <- duplicated(y$name.bi)
                    y[!dup, ]})})
uu <- lapply(uu, function(x) {
             rms <- lapply(x, nrow) == 0
             x[!rms]})
uu <- lapply(uu, function(x) Reduce(rbind, x))
dat <- Reduce(rbind, uu)
dat <- dat[-grep('[0-9\\-\\.]', dat$name.bi, perl = TRUE), ]
