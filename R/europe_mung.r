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

# no missing temporal informaiton
#
# formation names
# problems
#   roman numerals
#   A, B, etc. versions
#   alphanumeric seperators

eur$NAME <- gsub(pattern = '[0-9. ]*$', replacement = '', 
                 x = eur$NAME, perl = TRUE)  # trailing numbers and spaces
eur$NAME <- gsub(pattern = '\\([^)]*\\)', replacement = '', 
                 x = eur$NAME, perl = TRUE)  # parens

# trailing and leading spaces
eur$NAME <- gsub('^\\s+', '', eur$NAME)
eur$NAME <- gsub('\\s+$', '', eur$NAME)


# exclude aquative and volant nonsense
aq <- c('Cetacea', 'Desmostylia', 'Sirenia', 'Chiroptera')
eur <- eur[!(eur$ORDER %in% aq), ]
