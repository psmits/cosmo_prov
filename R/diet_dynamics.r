library(igraph)
library(plyr)
library(parallel)

source('../R/na_mung.r')
source('../R/europe_mung.r')

source('../R/bin_network.r')
source('../R/biogeo_struct.r')
source('../R/biogeo_bootstrap.r')
source('../R/window.r')

diet <- split(dat, f = dat$comdiet)
eurdt <- split(eur, f = eur$DIET_2)

stdiet <- lapply(diet, function(x) {
                 split(x, x$stage)})

stdigr <- lapply(stdiet, function(x) {
                 lapply(x, bin.network, taxa = 'name.bi', loc = 'formation')})
stdigr.bg <- lapply(stdigr, function(x) {
                    lapply(biogeosum, function(y) lapply(x, y))})

stdigr.hier <- lapply(stdigr, function(x) {
                      lapply(x, get.hier, level = 'family_name', data = dat)})
#stdigr.boot <- Map(function(x, y) {
#                   mclapply(biogeosum, function(foo) {
#                            mapply(biogeo.boot,
#                                   graph = x, taxon = y,
#                                   MoreArgs = list(fun = foo, 
#                                                   data = dat, 
#                                                   nsim = 10),
#                                   SIMPLIFY = FALSE)}, 
#                            mc.cores = detectCores())},
#                   x = stdigr, y = stdigr.hier)

# with explicit bins
wdth <- 2
dietwin <- lapply(diet, function(x) {
                  network.bin(x, width = wdth, time = 'ma_mid',
                              taxa = 'name.bi', loc = 'formation')})
dtwin.bg <- lapply(dietwin, function(x) {
                   lapply(biogeosum, function(y) lapply(x, y))})

dtwin.hier <- lapply(dietwin, function(x) {
                     lapply(x, get.hier, level = 'family_name', data = dat)})

dteur <- lapply(eurdt, function(x) {
                  network.bin(x, width = wdth, time = 'MID_AGE',
                              taxa = 'name.bi', loc = 'NAME')})
dteur.bg <- lapply(dteur, function(x) {
                   lapply(biogeosum, function(y) lapply(x, y))})
#dtwin.boot <- Map(function(x, y) {
#                   mclapply(biogeosum, function(foo) {
#                            mapply(biogeo.boot,
#                                   graph = x, taxon = y,
#                                   MoreArgs = list(fun = foo, 
#                                                   data = dat, 
#                                                   nsim = 10),
#                                   SIMPLIFY = FALSE)},
#                            mc.cores = detectCores())},
#                   x = dietwin, y = dtwin.hier)
