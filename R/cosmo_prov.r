library(igraph)
library(plyr)
library(parallel)

source('../R/na_mung.r')
source('../R/europe_mung.r')

source('../R/biogeo_struct.r')
source('../R/bin_network.r')
source('../R/window.r')

biogeo <- function(taxawin) {
  biocom <- lapply(taxawin, infomap.community)
  biomes <- Map(contract.vertices, taxawin, lapply(biocom, membership))
  biomes <- lapply(biomes, function(x){
                   x$weight = 1
                   x})
  biomes <- lapply(biomes, simplify)
  mem <- lapply(biocom, function(x) x$membership)

  win.bg <- list(bc = lapply(taxawin, bc),
                 end = Map(endemic, graph = taxawin, membership = mem),
                 avgcoc = Map(avgocc, graph = taxawin, membership = mem),
                 code = lapply(taxawin, code))
  win.bg
}

taxawin <- network.bin(dat, bin = 'bins', taxa = 'name.bi', loc = 'gid')
win.bg <- biogeo(taxawin)  # north america

eurwin <- network.bin(eur, bin = 'bins', taxa = 'name.bi', loc = 'gid')
eurwin.bg <- biogeo(eurwin)  # europe


# generic level
naocc.mat <- unique(dat[, c('occurrence.genus_name', 'bins', 'gid')])
nagenwin <- network.bin(naocc.mat, bin = 'bins', 
                        taxa = 'occurrence.genus_name', loc = 'gid')
erocc.mat <- unique(eur[, c('occurrence.genus_name', 'bins', 'gid')])
ergenwin <- network.bin(erocc.mat, bin = 'bins', 
                        taxa = 'occurrence.genus_name', loc = 'gid')


# need to do this for traits
biogeo.trait <- function(taxawin, dat) {
  biocom <- lapply(taxawin, infomap.community)
  biomes <- Map(contract.vertices, taxawin, lapply(biocom, membership))
  biomes <- lapply(biomes, function(x){
                   x$weight = 1
                   x})
  biomes <- lapply(biomes, simplify)
  mem <- lapply(biocom, function(x) x$membership)

  diet <- split(dat, f = dat$comdiet)
  dietwin <- lapply(diet, function(x) {
                    network.bin(x, bin = 'bins', taxa = 'name.bi', loc = 'gid')})
  life <- split(dat, f = dat$comlife)
  lifewin <- lapply(life, function(x) {
                    network.bin(x, bin = 'bins', taxa = 'name.bi', loc = 'gid')})

  trait <- lapply(taxawin, function(x) {
                  mm <- dat$name.bi %in% V(x)$name
                  unique(dat[mm, c('name.bi', 'comdiet', 'comlife')])})
  trait <- lapply(trait, function(x) x[!(duplicated(x$name.bi)), ])
  dtwin.bg <- list(bc = lapply(dietwin, function(x) lapply(x, bc)),
                   avgcoc = Map(avgocc, graph = taxawin, 
                                membership = mem, 
                                trait = lapply(trait, function(x) x$comdiet)),
                   end = Map(endemic, graph = taxawin, 
                             membership = mem, 
                             trait = lapply(trait, function(x) x$comdiet)),
                   code = lapply(dietwin, function(x) lapply(x, code)))
  lfwin.bg <- list(bc = lapply(lifewin, function(x) lapply(x, bc)),
                   avgcoc = Map(avgocc, graph = taxawin, 
                                membership = mem, 
                                trait = lapply(trait, function(x) x$comlife)),
                   end = Map(endemic, graph = taxawin, 
                             membership = mem, 
                             trait = lapply(trait, function(x) x$comlife)),
                   code = lapply(lifewin, function(x) lapply(x, code)))

  out <- list(diet = dtwin.bg, life= lfwin.bg)
  out
}

na.trait <- biogeo.trait(taxawin, dat)  # north america
er.trait <- biogeo.trait(eurwin, eur)  # europe
