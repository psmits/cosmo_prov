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
nagen.bg <- biogeo(nagenwin) # generic na

erocc.mat <- unique(eur[, c('occurrence.genus_name', 'bins', 'gid')])
ergenwin <- network.bin(erocc.mat, bin = 'bins', 
                        taxa = 'occurrence.genus_name', loc = 'gid')
ergen.bg <- biogeo(ergenwin) # generic europe

save.image(file = '../data/cosmo_prov_setup.r')
