library(igraph)
library(plyr)
library(parallel)

source('../R/na_mung.r')
source('../R/europe_mung.r')

source('../R/bin_network.r')
source('../R/biogeo_struct.r')
source('../R/window.r')

diet <- split(dat, f = dat$comdiet)
eurdt <- split(eur, f = eur$comdiet)

# with explicit bins
wdth <- 2
dietwin <- lapply(diet, function(x) {
                  network.bin(x, width = wdth, time = 'ma_mid',
                              taxa = 'name.bi', loc = 'formation')})
dtwin.bg <- lapply(dietwin, function(x) {
                   lapply(biogeosum, function(y) lapply(x, y))})

dteur <- lapply(eurdt, function(x) {
                  network.bin(x, width = wdth, time = 'ma_mid',
                              taxa = 'name.bi', loc = 'formation')})
dteur.bg <- lapply(dteur, function(x) {
                   lapply(biogeosum, function(y) lapply(x, y))})
