library(plyr)

source('../R/cosmo_prov.r')
source('../R/oxygen_curve.r')

# autocorrelation
na.acf <- lapply(win.bg, function(x) acf(unlist(x), 
                                         plot = FALSE, 
                                         na.action = na.pass))
er.acf <- lapply(eurwin.bg, function(x) acf(unlist(x), 
                                            plot = FALSE, 
                                            na.action = na.pass))

# correlations with the bins between the four different biogeo stats and oxygen
# using first differences
naoxym <- oo[names(oo) %in% names(taxawin)]
na.cor <- lapply(win.bg, function(x) cor.test(diff(rev(unlist(x))),
                                              diff(rev(naoxym)), 
                                              method = 'spearman',
                                              exact = FALSE))
eroxym <- oo[names(oo) %in% names(eurwin)]
er.cor <- lapply(eurwin.bg, function(x) cor.test(diff(rev(unlist(x))),
                                                 diff(rev(eroxym)), 
                                                 method = 'spearman',
                                                 exact = FALSE))

# correlations between continents
# using first differences
same <- Map(function(x, y) names(x) %in% names(y), x = win.bg, eurwin.bg)
nacomp <- Map(function(x, y) x[y], win.bg, same)
reg.cor <- Map(function(x, y) cor.test(diff(rev(unlist(x))),
                                       diff(rev(unlist(y))),
                                       method = 'spearman',
                                       exact = FALSE), 
               x = nacomp, y = eurwin.bg)
