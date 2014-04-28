library(igraph)
library(plyr)
library(parallel)

source('../R/na_mung.r')
source('../R/europe_mung.r')

source('../R/biogeo_struct.r')
source('../R/bin_network.r')
source('../R/window.r')
#source('../R/phylo_gen.r')

# with explicit bins
wdth <- 2
taxawin <- network.bin(dat, bin = 'bins', taxa = 'name.bi', loc = 'gid')
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

eurwin <- network.bin(eur, bin = 'bins', taxa = 'name.bi', loc = 'gid')
eurcom <- lapply(eurwin, infomap.community)
eurmes <- Map(contract.vertices, eurwin, lapply(eurcom, membership))
eurmes <- lapply(eurmes, function(x){
                 x$weight = 1
                 x})
eurmes <- lapply(eurmes, simplify)
mem <- lapply(eurcom, function(x) x$membership)
eur.bg <- list(bc = lapply(eurwin, bc),
               end = Map(endemic, graph = eurwin, membership = mem),
               avgcoc = Map(avgocc, graph = eurwin, membership = mem),
               code = lapply(eurwin, code))
