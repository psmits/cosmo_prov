## insert copy statement here

require(ggplot2)
require(scales)
require(mgcv)
require(reshape2)

source('../R/mammal_curve.R')

## generate zachos et al. 2008 curve
zachos <- read.csv('../data/2008_zachos_data.csv', header = T)
tits <- c('site', 'age', 'genus', 'd18o', 'd13c', 'd18o.5', 'd13c.5')
zachos <- zachos[, 1:7]
names(zachos) <- tits

zac.curve <- ggplot(zachos, aes(x = age, y = d18o)) + geom_point(aes(alpha = 0.5))
zac.curve <- zac.curve + stat_smooth(method = 'gam', 
                                     formula = y ~ s(x, bs = 'tp'), 
                                     na.rm = T)
zac.curve <- zac.curve + scale_y_reverse()
zac.curve <- zac.curve + theme(panel.grid = element_blank())
zac.curve <- zac.curve + theme(panel.background = element_blank())
zac.curve <- zac.curve + theme(panel.border = element_blank())
zac.curve <- zac.curve + theme(axis.line = element_line(colour = 'black'))

zac.gam <- gam(d18o ~ s(age, bs = 'tp'), data = zachos, method = 'REML')

## make a binned version of the curve at the million year level
obins <- 1:66
ass.oxy <- vector(mode = 'list', length = (length(obins) - 1))
for (kk in seq(length(ass.oxy))) {
  ass.oxy[[kk]] <- which(zachos$age >= obins[kk] & zachos$age < obins[kk + 1])
}

avger <- function(pos, z) {
  mean(z$d18o[pos], na.rm = TRUE)
}

oxy.bins <- melt(lapply(ass.oxy, avger, zachos))
colnames(oxy.bins) <- c('d18', 'time')

mid.ox <- c()
for (ii in seq(nrow(bins))) {
  ones <- zachos$age > bins.2[ii, 2] & zachos$age < bins.2[ii, 1]
  mid.ox[ii] <- mean(zachos$d18o[ones], na.rm = T)
}

mid.bin <- mammal2s$mid
zac.mam.bin <- data.frame(mid.bin, mid.ox)

zac.bc <- ggplot(zac.mam.bin, aes(x = mid.bin, y = mid.ox)) + geom_point()
zac.bc <- zac.bc + geom_point() + stat_smooth(na.rm = T)
zac.bc <- zac.bc + scale_y_reverse()
zac.bc <- zac.bc + scale_x_continuous(limits = c(0, 70))

