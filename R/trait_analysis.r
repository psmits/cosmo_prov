library(plyr)
library(reshape2)

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
# there are more locomotor categories than diet categories
io <- list()
for(ii in seq(length(dtwin.bg))) {
  dit <- dtwin.bg[[ii]]
  out <- list()
#  lapply(dit, function(x) print(length(x)))
  for(jj in seq(length(lfwin.bg))) {
    ll <- lfwin.bg[[jj]]
#    lapply(ll, function(x) print(length(x)))
    out[[jj]] <- mapply(function(x, y) {
                        x <- unlist(x)
                        y <- unlist(y)
                        xx <- x[names(x) %in% names(y)]
                        yy <- y[names(y) %in% names(x)]
                        try(cor.test(diff(rev(xx)), diff(rev(yy)), 
                                     method = 'spearman', 
                                     exact = FALSE))
                   }, x = dit, y = ll, SIMPLIFY = FALSE)
  }
  names(out) <- names(lfwin.bg)
  io[[ii]] <- out
}
names(io) <- names(dtwin.bg)
