library(plyr)
library(paleoTS)

source('../R/cosmo_prov.r')
source('../R/oxygen_curve.r')

# correlations with the bins between the four different biogeo stats and oxygen
oxym <- oo[names(oo) %in% names(taxawin)]
win.cor <- lapply(win.bg, function(x) cor.test(diff(rev(unlist(x))),
                                               diff(rev(oxym)), 
                                               method = 'spearman',
                                               exact = FALSE))

# paleoTS using the bootstrap results from each bin
tx.bt <- lapply(taxawin.boot, function(x) {
                lapply(x, function(y) {
                       list(mm = mean(y), 
                            vv = var(y),
                            nn = length(y))})})
tt <- as.numeric(names(taxawin.boot$bc))
mtb <- melt(tx.bt)       
mtb <- split(mtb, mtb$L3)
mtb <- lapply(mtb, function(x) split(x, x$L1))
mtb.ts <- mapply(function(x, y, z) cbind(mm = x$value,
                                         vv = y$value,
                                         nn = z$value,
                                         tt = tt), 
                 x = mtb$mm, y = mtb$vv, z = mtb$nn,
                 SIMPLIFY = FALSE)
mtb.ts <- lapply(mtb.ts, as.data.frame)
mtb.pts <- lapply(mtb.ts, function(x) as.paleoTS(mm = x$mm, 
                                                 vv = x$vv, 
                                                 nn = x$nn, 
                                                 tt = x$tt))
mtb.res <- lapply(mtb.pts, fit3models, silent = TRUE, method = 'Joint')
