library(plyr)

source('../R/diet_dynamics.r')
source('../R/life_dynamics.r')
source('../R/oxygen_curve.r')

# diet ~ oxygen
oxyd <- Map(function(x, y) {
            x[names(x) %in% names(y)]
            x}, x = oo, y = dietwin)
oxyd <- unlist(oxyd)
diet.cor <- lapply(dtwin.bg, function(x) {
                   out <- list()
                   for (ii in seq(length(x))) {
                     val <- unlist(x[[ii]])
                     ox <- oxyd[names(oxyd) %in% names(val)]
                     out[[ii]] <- cor.test(diff(rev(val)), diff(rev(ox)), 
                                           method = 'spearman',
                                           exact = FALSE)
                   }
                   out})



# locomotor ~ oxygen
loco.cor <- lapply(lfwin.bg, function(x) {
                   out <- list()
                   for (ii in seq(length(x))) {
                     val <- unlist(x[[ii]])
                     ox <- oxyd[names(oxyd) %in% names(val)]
                     out[[ii]] <- try(cor.test(diff(rev(val)), diff(rev(ox)), 
                                               method = 'spearman',
                                               exact = FALSE))
                   }
                   names(out) <- names(x)
                   out})

# diet ~ locomotor
# compare like with like
