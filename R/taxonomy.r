library(taxize)
library(plyr)
library(stringr)

source('../R/my_get_eolid.r')
source('../R/eol_taxonomy.r')
source('../R/taxon_names.r')

source('../R/extinct_families.r')
source('../R/na_mung.r')
source('../R/europe_mung.r')

eol.key = '2a9932f264f3f0421db36158b6e785b535c6da0e'
now <- read.csv('../data/now_database.csv')
tomiya <- read.csv('../data/tomiya_mass.csv', stringsAsFactors = FALSE)

# update from the now
miss.fam <- dat$occurrence.genus_name[dat$family_name == '']
in.now <- unique(miss.fam[miss.fam %in% now$GENUS])
now.info <- now[now$GENUS %in% in.now, ]
now.info <- now.info[now.info$FAMILY != 'incertae sedis', ]
now.heir <- unique(now.info[, c('ORDER', 'FAMILY', 'GENUS')])
dat <- replace.taxonomy(dat, now.heir)

miss.fam <- eur$occurrence.genus_name[eur$family_name == '']
in.now <- unique(miss.fam[miss.fam %in% now$GENUS])
now.info <- now[now$GENUS %in% in.now, ]
now.info <- now.info[now.info$FAMILY != 'incertae sedis', ]
now.heir <- unique(now.info[, c('ORDER', 'FAMILY', 'GENUS')])
eur <- replace.taxonomy(eur, now.heir)

# update from Susumu
miss.fam <- dat$occurrence.genus_name[dat$family_name == '']
in.tom <- unique(miss.fam[miss.fam %in% tomiya$Genus])
tom.info <- tomiya[tomiya$Genus %in% in.tom, ]
tom.info <- tom.info[tom.info$Family != 'unknown', ]
tom.heir <- unique(tom.info[, c('Order', 'Family', 'Genus')])
dat <- replace.taxonomy(dat, tom.heir)

miss.fam <- eur$occurrence.genus_name[eur$family_name == '']
in.tom <- unique(miss.fam[miss.fam %in% tomiya$Genus])
tom.info <- tomiya[tomiya$Genus %in% in.tom, ]
tom.info <- tom.info[tom.info$Family != 'unknown', ]
tom.heir <- unique(tom.info[, c('Order', 'Family', 'Genus')])
eur <- replace.taxonomy(eur, tom.heir)


# update with EOL
genera <- unique(c(dat$occurrence.genus_name, eur$occurrence.genus_name))
genera <- genera[genera != 'Vulpes']
heir <- grab.heir(genera, key = eol.key)  # updated hierarchy information
orders <- col_downstream(name = 'Mammalia', downto = 'Order')
new.orders <- c('Eulipotyphla', # combines erinaceidae and soricomorpha
                'Cetartiodactyla' # combines cetacea and artiodactyla
                )
orders <- c(orders$Mammalia[, 2], new.orders)
heir <- heir[(heir[, 1] %in% orders), ]  # there are random plants/insects
# for matching genera, update family and order information
udat <- replace.taxonomy(dat, heir)
ueur <- replace.taxonomy(eur, heir)
udat <- miss.ord(udat)
ueur <- miss.ord(ueur)


# extinct families
extinct.heir <- Reduce(rbind, Map(function(x, y) cbind(order = rep(x, length(y)), 
                                                       family = y),
                                  x = names(extinct), y = extinct))
bonus.ord <- unique(c(orders, names(extinct)))

replace.extinct <- function(data, heir) {
  for(ii in seq(nrow(heir))) {
    rr <- which(data$family_name == heir[ii, 2])
    data[rr, 'order_name'] <- heir[ii, 1]
  }
  data$order_name[data$order_name == 'Acreodi'] <- 'Mesonychia'
  data
}
udat <- replace.extinct(udat, extinct.heir)
ueur <- replace.extinct(ueur, extinct.heir)
udat$family_name[grep('[0-9]', udat$family_name)] <- ''
udat$order_name[grep('[0-9]', udat$order_name)] <- ''
ueur$family_name[grep('[0-9]', ueur$family_name)] <- ''
ueur$order_name[grep('[0-9]', ueur$order_name)] <- ''

keep <- c('order_name', 'family_name', 'occurrence.genus_name', 'name.bi', 
          'ma_min', 'ma_max')
na.tax <- udat[, keep]
er.tax <- ueur[, keep]

save(na.tax, er.tax, 
     file = '../data/update_taxonomy.rdata')
