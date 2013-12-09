library(igraph)
library(plyr)
library(parallel)

source('../R/na_mung.r')
source('../R/europe_mung.r')

source('../R/bin_network.r')
source('../R/biogeo_struct.r')
source('../R/biogeo_bootstrap.r')
source('../R/window.r')

stage <- split(dat, f = dat$stage)
stgraph <- lapply(stage, bin.network, taxa = 'name.bi', loc = 'formation')

stgraph.bg <- lapply(biogeosum, function(x) {
                     lapply(stgraph, x)})

#stgraph.hier <- lapply(stgraph, get.hier, level = 'family_name', data = dat)
#stgraph.boot <- mclapply(biogeosum, function(foo) {
#                         mapply(biogeo.boot,
#                         graph = stgraph, taxon = stgraph.hier,
#                         MoreArgs = list(fun = foo, data = dat, nsim = 10),
#                         SIMPLIFY = FALSE)}, mc.cores = detectCores())

# with explicit bins
wdth <- 2
taxawin <- network.bin(dat, width = wdth, time = 'ma_mid', 
                         taxa = 'name.bi', loc = 'formation')
win.bg <- lapply(biogeosum, function(x) {
                 lapply(taxawin, x)})

eurwin <- network.bin(eur, width = wdth, time = 'ma_mid',
                      taxa = 'name.bi', loc = 'formation')

eurwin.bg <- lapply(biogeosum, function(x) {
                    lapply(eurwin, x)})

#taxawin.hier <- lapply(taxawin, get.hier, level = 'family_name', data = dat)
#taxawin.boot <- mclapply(biogeosum, function(foo) {
#                         mapply(biogeo.boot,
#                         graph = taxawin, taxon = taxawin.hier,
#                         MoreArgs = list(fun = foo, data = dat, nsim = 10),
#                         SIMPLIFY = FALSE)}, mc.cores = detectCores())
