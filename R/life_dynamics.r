library(igraph)
library(plyr)
library(parallel)

#source('../R/na_mung.r')
#source('../R/europe_mung.r')

source('../R/bin_network.r')
source('../R/biogeo_struct.r')
source('../R/window.r')

life <- split(dat, f = dat$comlife)
eulf <- split(eur, f = eur$comlife)

# with explicit bins
wlfh <- 2
lifewin <- lapply(life, function(x) {
                  network.bin(x, bin = 'bins', taxa = 'name.bi', loc = 'gid')})
lfwin.bg <- lapply(lifewin, function(x) {
                   lapply(biogeosum, function(y) lapply(x, y))})

lfeur <- lapply(eulf, function(x) {
                network.bin(x, bin = 'bins', taxa = 'name.bi', loc = 'gid')})
lfeur.bg <- lapply(lfeur, function(x) {
                   lapply(biogeosum, function(y) lapply(x, y))})
