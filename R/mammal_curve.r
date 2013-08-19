## make the diversity curves for mammal2ss
require(ggplot2)
require(scales)
require(reshape2)

prop.name <- c('Bin', 'name', 'base', 'mid', 'n.collec', 'mean.div', 'origination', 'extinction')

# 2 million year diversity
mammal2s <- read.csv('../data/psmamall_subsample.csv', header = T)
names(mammal2s) <- prop.name

bin.widths <- mammal2s$base - mammal2s$mid
bots <- mammal2s$base - 2 * bin.widths
bins <- cbind(mammal2s$base, bots)

colnames(bins) <- c('old', 'young')
bins.2 <- as.data.frame(bins)  # useful value on export for the oxygen curves

# 10 million year diversity 
#mammal10s <- read.csv('psmamocc_10my_subsampled_curve_data.csv', header = T)
#names(mammal10s) <- prop.name

# subepoch diversity
#mammalsubs <- read.csv('psmamocc_subep_subsampled_curve_data.csv', header = T)
#names(mammalsubs) <- prop.name

gmam2 <- ggplot(mammal2s, aes(x = mid, y = mean.div)) + geom_line() + stat_smooth()
gmam2 <- gmam2 + scale_x_continuous(name = 'Time (Mya)')
gmam2 <- gmam2 + scale_y_continuous(name = 'Generic diversity'
                                    )

gmam2o <- ggplot(mammal2s, aes(x = mid, y = origination))
gmam2o <- gmam2o + geom_line() + stat_smooth()
gmam2o <- gmam2o + scale_x_continuous(name = 'Time (Mya)'
                                      , limits = c(0, 70)
                                      )
gmam2o <- gmam2o + scale_y_continuous(name = 'Origination rate'
                                      , trans = log10_trans()
                                      )

gmam2e <- ggplot(mammal2s, aes(x = mid, y = extinction)) + geom_line() + stat_smooth()
gmam2e <- gmam2e + scale_x_continuous(name = 'Time (Mya)'
                                      , limits = c(0, 70)
                                      )
gmam2e <- gmam2e + scale_y_continuous(name = 'Extinction rate'
                                      , trans = log10_trans()
                                      )

