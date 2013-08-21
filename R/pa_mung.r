library(plyr)
library(reshape2)

source('../R/clean_pbdb.r')
source('../R/collapse_names.r')
source('../R/remove_zeroes.r')
source('../R/taxon_names.r')

dat <- read.csv('../data/mam-occs.csv', stringsAsFactors = FALSE)

# if missing any temporal information, exclude
temp <- c('period', 'epoch', 'subepoch', 'stage', 'formation')
trm <- Reduce(union, lapply(temp, function(x) which(dat[, x] == '')))
dat <- dat[-trm, ]

# clean up all the formation names
rmform <- c('middle Miocene', 'troublesome', 'vertebrate', 'unnamed sandstone')
dat <- dat[!dat$formation %in% rmform, ]
dat$formation  <- gsub(pattern = '[\\"?]', 
                       replacement = '', 
                       x = dat$formation, perl = TRUE)
dat$formation <- gsub(pattern = '\\([^)]*\\)', 
                      replacement = '', 
                      x = dat$formation, 
                      perl = TRUE)
dat$formation <- gsub('^[a-z. ]*', '', dat$formation, perl = TRUE)
dat$formation <- sub('^\\s+', '', dat$formation)
dat$formation <- sub('\\s+$', '', dat$formation)
dat$formation[grep('Fort Union', dat$formation)] <- 'Fort Union'

# make the number columns numeric
num <- c('max_ma', 'max_ma_error',
         'min_ma', 'min_ma_error',
         'interval_base', 'interval_top', 'interval_midpoint',
         'interpolated_base', 'interpolated_top', 'interpolated_mid',
         'ma_max', 'ma_min', 'ma_mid',
         'collection_no')
dat[, num] <- apply(dat[, num], 2, as.numeric)
# need to have ma_mid values
dat <- dat[!is.na(dat$ma_mid), ]

# make sure i'm entirely in the Cenozoic via max
kpg <- 65.6
if(sum(dat$interval_top > kpg, na.rm = TRUE)) {
  dat <- dat[!dat$interval_top > kpg]
}

# exclude aquatic and volant nonsense
aq <- c('Cetacea', 'Desmostylia', 'Sirenia', 'Chiroptera')
dat <- dat[!(dat$order_name %in% aq), ]

lf <- c('amphibious', 'volant', 'aquatic')
dat <- dat[!(dat$life_habit %in% lf), ]

# remove all the sp.-s
#grep('sp', x = dat$occurence.species_name, perl = TRUE)
dat <- dat[dat$occurrence.species_name != 'sp.', ]  # change to a grep 


binm <- with(dat, binom.make(occurrence.genus_name, occurrence.species_name))
dat <- cbind(dat, name.bi = binm)
dat$name.bi <- as.character(dat$name.bi)


# there are a lot of dupilcated taxa in each formation
# split data set into formation and then remove duplicate taxa

# but first, grab the complete version
pa.mat <- dat
dd <- split(dat, dat$formation)
dd <- lapply(dd, function(x) {
             x <- x[!duplicated(x$name.bi), ]
             x})
dat <- Reduce(rbind, dd)

# make a rougher diet column 
herb <- c('herbivore', 'grazer', 'browser', 'folivore', 'granivore')
omm <- c('frugivore', 'omnivore')
car <- c('carnivore', 'insectivore')
dat$comdiet <- dat$diet1
dat$comdiet[dat$diet1 %in% herb] <- 'herb'
dat$comdiet[dat$diet1 %in% omm] <- 'omni'
dat$comdiet[dat$diet1 %in% car] <- 'carni'
pa.mat$comdiet <- pa.mat$diet1
pa.mat$comdiet[pa.mat$diet1 %in% herb] <- 'herb'
pa.mat$comdiet[pa.mat$diet1 %in% omm] <- 'omni'
pa.mat$comdiet[pa.mat$diet1 %in% car] <- 'carni'


st <- split(dat, dat$stage)
yst <- lapply(st, function(x) mean(x$ma_mid))
yst <- Map(function(x, y) rep(x, y), x = yst, y = lapply(st, nrow))
st <- Map(function(x, y) cbind(y, stmid = x), x = yst, y = st)
dat <- Reduce(rbind, st)
