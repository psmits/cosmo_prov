library(igraph)
library(plyr)
library(parallel)

source('../R/bin_network.r')
source('../R/biogeo_struct.r')
source('../R/window.r')

diet <- split(dat, f = dat$comdiet)
eurdt <- split(eur, f = eur$comdiet)

# with explicit bins
wdth <- 2
dietwin <- lapply(diet, function(x) {
                  network.bin(x, bin = 'bins', taxa = 'name.bi', loc = 'gid')})
dtwin.bg <- lapply(dietwin, function(x) {
                   lapply(biogeosum, function(y) lapply(x, y))})

dteur <- lapply(eurdt, function(x) {
                  network.bin(x, bin = 'bins', taxa = 'name.bi', loc = 'gid')})
dteur.bg <- lapply(dteur, function(x) {
                   lapply(biogeosum, function(y) lapply(x, y))})
