library(taxize)
library(plyr)

source('../R/my_get_eolid.r')
source('../R/eol_taxonomy.r')
source('../R/taxon_names.r')
source('../R/phylo_gen.r')

source('../R/extinct_families.r')
source('../R/na_mung.r')
source('../R/europe_mung.r')

eol.key = '2a9932f264f3f0421db36158b6e785b535c6da0e'

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
  for(ii in seq(nrow(extinct.heir))) {
    rr <- which(data$family_name == heir[ii, 2])
    data[rr, 'order_name'] <- heir[ii, 1]
  }
  data$order_name[data$order_name == 'Acreodi'] <- 'Mesonychia'
  data
}
udat <- replace.extinct(udat, extinct.heir)
ueur <- replace.extinct(ueur, extinct.heir)

no.fam <- which(udat$family_name == '')
no.ord <- which(udat$order_name == '')
rms <- unique(c(no.fam, no.ord))
clean.na <- udat[-rms, ]
na.tree <- big.tree(clean.na[, c('occurrence.genus_name', 
                                 'family_name', 
                                 'order_name', 
                                 'name.bi')])

no.fam <- which(ueur$family_name == '')
no.ord <- which(ueur$order_name == '')
rms <- unique(c(no.fam, no.ord))
clean.er <- ueur[-rms, ]
er.tree <- big.tree(clean.er[, c('occurrence.genus_name', 
                                 'family_name', 
                                 'order_name', 
                                 'name.bi')])
