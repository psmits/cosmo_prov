library(igraph)
library(plyr)
library(parallel)

source('../R/na_mung.r')

source('../R/bin_network.r')
source('../R/biogeo_struct.r')
source('../R/biogeo_bootstrap.r')
source('../R/window.r')

life <- split(dat, f = dat$life_habit)

stlife <- lapply(life, function(x) {
                 split(x, x$stage)})

stlfgr <- lapply(stlife, function(x) {
                 lapply(x, bin.network, taxa = 'name.bi', loc = 'formation')})
stlfgr.bg <- lapply(stlfgr, function(x) {
                    lapply(biogeosum, function(y) lapply(x, y))})

stlfgr.hier <- lapply(stlfgr, function(x) {
                      lapply(x, get.hier, level = 'family_name', data = dat)})
stlfgr.boot <- Map(function(x, y) {
                   mclapply(biogeosum, function(foo) {
                            mapply(biogeo.boot,
                                   graph = x, taxon = y,
                                   MoreArgs = list(fun = foo, 
                                                   data = dat, 
                                                   nsim = 10),
                                   SIMPLIFY = FALSE)},
                            mc.cores = detectCores())},
                   x = stlfgr, y = stlfgr.hier)

# with explicit bins
wlfh <- 2
lifewin <- lapply(life, function(x) {
                  network.bin(x, width = wlfh, time = 'ma_mid',
                              taxa = 'name.bi', loc = 'formation')})
lfwin.bg <- lapply(lifewin, function(x) {
                   lapply(biogeosum, function(y) lapply(x, y))})

lfwin.hier <- lapply(lifewin, function(x) {
                     lapply(x, get.hier, level = 'family_name', data = dat)})
lfwin.boot <- Map(function(x, y) {
                   mclapply(biogeosum, function(foo) {
                            mapply(biogeo.boot,
                                   graph = x, taxon = y,
                                   MoreArgs = list(fun = foo, 
                                                   data = dat, 
                                                   nsim = 10),
                                   SIMPLIFY = FALSE)},
                            mc.cores = detectCores())},
                   x = lifewin, y = lfwin.hier)
