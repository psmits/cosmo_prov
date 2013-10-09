library(plyr)
library(reshape2)

source('../R/clean_pbdb.r')
source('../R/collapse_names.r')
source('../R/remove_zeroes.r')
source('../R/taxon_names.r')

eur <- read.csv('../data/euro-occs.csv', stringsAsFactors = FALSE)

# remove things with missing temporal information
temp <- c('period', 'epoch', 'subepoch', 'stage', 'formation')
trm <- Reduce(union, lapply(temp, function(x) which(eur[, x] == '')))
eur <- eur[-trm, ]

# clean up the formation names
# doesn't start with a capital letter
eur <- eur[grep(pattern = '^[A-Z]', eur$formation, perl = TRUE), ]
eur$formation <- gsub(pattern = ' s', replacement = ' S',
                      eur$formation, perl = TRUE)

# make the number columns numeric
num <- c('max_ma', 'max_ma_error',
         'min_ma', 'min_ma_error',
         'interval_base', 'interval_top', 'interval_midpoint',
         'interpolated_base', 'interpolated_top', 'interpolated_mid',
         'ma_max', 'ma_min', 'ma_mid',
         'collection_no')
eur[, num] <- apply(eur[, num], 2, as.numeric)
# need to have ma_mid values
eur <- eur[!is.na(eur$ma_mid), ]

# make sure i'm entirely in the Cenozoic
kpg <- 65.6
if(sum(eur$interval_top > kpg, na.rm = TRUE)) {
  eur <- eur[!eur$interval_top > kpg]
}

# aquatics and volants
aq <- c('Cetacea', 'Desmostylia', 'Sirenia', 'Chiroptera')
eur <- eur[!(eur$order_name %in% aq), ]

lf <- c('amphibious', 'volant', 'aquatic')
eur <- eur[!(eur$life_habit %in% lf), ]

# remove all the sp.-s
#grep('sp', x = dat$occurence.species_name, perl = TRUE)
eur <- eur[eur$occurrence.species_name != 'sp.', ]  # change to a grep 








#eur <- read.csv('../data/now_database.csv', stringsAsFactors = FALSE)

# restrict to europe
#eur.count <- c('Switzerland', 'Spain', 'Greece', 'Germany', 'Italy',
#               'France', 'Bulgaria', 'Ukraine', 'Australia', 'Portugal',
#               'Belgium', 'Romania', 'Moldova', 'United Kingdom', 'Hungary',
#               'Czech Republic', 'Poland', 'Serbia', 'Slovakia', 'Netherlands',
#               'Ireland', 'Croatia', 'Malta', 'Serbia and Montenegro', 'Sweden',
#               'Belarus', 'Slovenia', 'Estonia')
#eur <- eur[eur$COUNTRY %in% eur.count, ]

# restrict to Cenozoic
#kpg <- 65.6
#eur <- eur[eur$MAX_AGE <= kpg, ]

# clean locality names
#eur$NAME <- gsub(pattern = '[0-9](.*)$', replacement = '',
#                 x = eur$NAME, perl = TRUE)
#eur$NAME <- gsub(pattern = '\\([^)]*\\)', replacement = '', 
#                 x = eur$NAME, perl = TRUE)  # parens
#eur$NAME <- gsub(pattern = '[\\[\\(](.*)$', replacement = '', 
#                 x = eur$NAME, perl = TRUE)
# trailing and leading spaces
#eur$NAME <- gsub('^\\s+', '', eur$NAME)
#eur$NAME <- gsub('\\s+$', '', eur$NAME)
#eur$NAME <- gsub(pattern = '[-A-Z]$', replacement = '',
#                 x = eur$NAME, perl = TRUE)
#eur$NAME <- gsub(pattern = ' I*$', replacement = '', 
#                 x = eur$NAME, perl = TRUE)
#eur$NAME <- gsub(pattern = ' cave(.)+$', replacement = '', 
#                 x = eur$NAME, perl = TRUE)
#eur$NAME <- gsub(pattern = ' st(.)+$', replacement = '', 
#                 x = eur$NAME, perl = TRUE)
#eur$NAME <- gsub('^\\s+', '', eur$NAME)
#eur$NAME <- gsub('\\s+$', '', eur$NAME)


# exclude aquative and volant nonsense
#aq <- c('Cetacea', 'Desmostylia', 'Sirenia', 'Chiroptera')
#eur <- eur[!(eur$ORDER %in% aq), ]

# get read of any locality with only one entry
# these are useless
#rms <- which(table(eur$NAME) == 1)
#eur <- eur[!(eur$NAME %in% names(rms)), ]

# add in midpoint age
#eur$MID_AGE <- (eur$MAX_AGE - eur$MIN_AGE) / 2

# remove indeterminate genera
#eur <- eur[eur$GENUS != 'indet.', ]

# add binomial name
#binm <- with(eur, binom.make(GENUS, SPECIES))
#eur <- cbind(eur, name.bi = binm)
#eur$name.bi <- as.character(eur$name.bi)

# get rid of the duplicates
#aa <- split(eur, eur$NAME)
#aa <- lapply(aa, function(x) {
#             x <- x[!duplicated(x$name.bi), ]
#             x})
#eur <- Reduce(rbind, aa)

# get rid of taxa with no dietary information
#dts <- c('omnivore', 'herbivore', 'carnivore')
#eur <- eur[eur$DIET_2 %in% dts, ]
