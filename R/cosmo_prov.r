library(igraph)
library(plyr)
library(parallel)

source('../R/na_mung.r')
source('../R/europe_mung.r')

source('../R/bin_network.r')
source('../R/biogeo_struct.r')
source('../R/window.r')
source('../R/phylo_gen.r')

# with explicit bins
wdth <- 2
taxawin <- network.bin(dat, bin = 'bins', taxa = 'name.bi', loc = 'gid')
win.bg <- lapply(biogeosum, function(x) {
                 lapply(taxawin, x)})
win.phy <- lapply(taxawin, mean.coph, data = dat)
win.bg$phy <- lapply(win.phy, mean)

eurwin <- network.bin(eur, bin = 'bins', taxa = 'name.bi', loc = 'gid')
eurwin.bg <- lapply(biogeosum, function(x) {
                    lapply(eurwin, x)})
eur.phy <- lapply(eurwin, mean.coph, data = eur)
eurwin.bg$phy <- lapply(eur.phy, mean)

