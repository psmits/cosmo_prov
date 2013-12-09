library(igraph)
library(plyr)
library(parallel)

source('../R/na_mung.r')
source('../R/europe_mung.r')

source('../R/bin_network.r')
source('../R/biogeo_struct.r')
source('../R/window.r')

# with explicit bins
wdth <- 2
taxawin <- network.bin(dat, bin = 'bins', taxa = 'name.bi', loc = 'gid')
win.bg <- lapply(biogeosum, function(x) {
                 lapply(taxawin, x)})

eurwin <- network.bin(eur, bin = 'bins', taxa = 'name.bi', loc = 'gid')
eurwin.bg <- lapply(biogeosum, function(x) {
                    lapply(eurwin, x)})

