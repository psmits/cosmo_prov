library(plyr)
library(reshape2)

source('../R/clean_pbdb.r')
source('../R/collapse_names.r')
source('../R/remove_zeroes.r')
source('../R/taxon_names.r')

eur <- read.csv('../data/now_database.csv', stringsAsFactors = FALSE)

# restrict to europe
eur.count <- c('Switzerland', 'Spain', 'Greece', 'Germany', 'Italy',
               'France', 'Bulgaria', 'Ukraine', 'Australia', 'Portugal',
               'Belgium', 'Romania', 'Moldova', 'United Kingdom', 'Hungary',
               'Czech Republic', 'Poland', 'Serbia', 'Slovakia', 'Netherlands',
               'Ireland', 'Croatia', 'Malta', 'Serbia and Montenegro', 'Sweden',
               'Belarus', 'Slovenia', 'Estonia')
eur <- eur[eur$COUNTRY %in% eur.count, ]

# restrict to Cenozoic
kpg <- 65.6
eur <- eur[eur$MAX_AGE <= kpg, ]

# clean locality names
eur$NAME <- gsub(pattern = '[0-9](.*)$', replacement = '',
                 x = eur$NAME, perl = TRUE)
eur$NAME <- gsub(pattern = '\\([^)]*\\)', replacement = '', 
                 x = eur$NAME, perl = TRUE)  # parens
eur$NAME <- gsub(pattern = '[\\[\\(](.*)$', replacement = '', 
                 x = eur$NAME, perl = TRUE)
# trailing and leading spaces
eur$NAME <- gsub('^\\s+', '', eur$NAME)
eur$NAME <- gsub('\\s+$', '', eur$NAME)
eur$NAME <- gsub(pattern = '[-A-Z]$', replacement = '',
                 x = eur$NAME, perl = TRUE)
eur$NAME <- gsub(pattern = ' I*$', replacement = '', 
                 x = eur$NAME, perl = TRUE)
eur$NAME <- gsub(pattern = ' cave(.)+$', replacement = '', 
                 x = eur$NAME, perl = TRUE)
eur$NAME <- gsub(pattern = ' st(.)+$', replacement = '', 
                 x = eur$NAME, perl = TRUE)
eur$NAME <- gsub('^\\s+', '', eur$NAME)
eur$NAME <- gsub('\\s+$', '', eur$NAME)


# exclude aquative and volant nonsense
aq <- c('Cetacea', 'Desmostylia', 'Sirenia', 'Chiroptera')
eur <- eur[!(eur$ORDER %in% aq), ]

# get read of any locality with only one entry
# these are useless
rms <- which(table(eur$NAME) == 1)
eur <- eur[!(eur$NAME %in% names(rms)), ]
