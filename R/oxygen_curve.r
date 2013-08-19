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
