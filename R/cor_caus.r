# analysis of the correlation and causation between biotic and abiotic factors 
# on cosmopolitan and provincial dynamics in Cenozoic mammals

library(reshape2)

source('../R/cosmo_prov.r')
source('../R/diet_dynamics.r')

wnbg <- melt(win.bg)
wnbg <- split(wnbg, f = stbg$L1)

dtwbg <- melt(dtwin.bg)
dtwbg <- split(dtwbg, f = dtwbg$L1)
dtwbg <- lapply(dtwbg, function(x) split(x, f = x$L2))
