library(mgcv)
library(plyr)
library(reshape2)
library(ggplot2)
library(scales)

source('../R/cosmo_prov.r')
source('../R/diet_dynamics.r')

theme_set(theme_bw())
cbp <- c('#E69F00', '#56B4E9', '#009E73', '#F0E442', 
         '#0072B2', '#D55E00', '#CC79A7')

time.match <- cbind(ma = dat$stmid, st = data.frame(dat$stage))
tm <- unique(time.match)

# total graph plots
# stage
st.boot <- melt(stgraph.boot)
st.boot <- cbind(st.boot, year = tm[match(st.boot$L2, tm[, 2]), 1])
gstb <- ggplot(st.boot, aes(x = year, y = value))
gstb <- gstb + geom_point()
gstb <- gstb + facet_wrap(~ L1, scales = 'free_y')
gstb <- gstb + stat_smooth()

# explicit bin
win.boot <- melt(taxawin.boot)
win.boot$L2 <- as.numeric(win.boot$L2)
gwinbt <- ggplot(win.boot, aes(x = L2, y = value))
gwinbt <- gwinbt + geom_point() 
gwinbt <- gwinbt + facet_wrap(~ L1, scales = 'free_y')
gwinbt <- gwinbt + stat_smooth()


# diet graphs
# stage
dtst.boot <- melt(stdigr.boot)
dtst.boot <- cbind(dtst.boot, year = tm[match(dtst.boot$L3, tm[, 2]), 1])
gdstb <- ggplot(dtst.boot, aes(x = year, y = value, colour = L1))
gdstb <- gdstb + geom_point()
gdstb <- gdstb + facet_wrap(~ L2, scales = 'free_y')
gdstb <- gdstb + stat_smooth(na.action = na.omit)

# explicit bin
dtw.bt <- melt(dtwin.boot)
dtw.bt$value[is.nan(dtw.bt$value)] <- NA
dtw.bt$value[is.infinite(dtw.bt$value)] <- NA
dtw.bt$L3 <- as.numeric(dtw.bt$L3)
gdwbt <- ggplot(dtw.bt, aes(x = L3, y = value, colour = L1))
gdwbt <- gdwbt + geom_point() 
gdwbt <- gdwbt + facet_wrap(~ L2, scales = 'free_y')
gdwbt <- gdwbt + stat_smooth(na.action = na.omit)
