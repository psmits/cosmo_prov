library(plyr)
library(reshape2)
library(mapproj)
library(stringr)
library(dismo)
library(raster)
library(sp)
library(XML)
library(maptools)
library(foreign)
library(rgdal)
library(ggplot2)
library(grid)

source('../R/na_mung.r')

theme_set(theme_bw())
cbp <- c('#E69F00', '#56B4E9', '#009E73', '#F0E442', 
         '#0072B2', '#D55E00', '#CC79A7')
theme_update(axis.text = element_text(size = 20),
             axis.title = element_text(size = 30),
             legend.text = element_text(size = 25),
             legend.title = element_text(size = 26),
             legend.key.size = unit(2, 'cm'),
             strip.text = element_text(size = 20))

data <- data.frame(id = rownames(globe.map@data),
                   values = sample(1:10, length(globe.map), replace = TRUE),
                   globe.map@data, stringsAsFactors = FALSE)

data_fort <- fortify(globe.map)
data_merged <- join(data_fort, data, by = 'id')
# make a map of the globe
# ggplot(data_merged, aes(x = long, y = lat, group = group)) + geom_path()

map.pts <- data.frame(rasterToPoints(sp.ras))
names(map.pts) <- c('long', 'lat', 'map')


# plot the occurrences
map.occ <- ggplot(map.pts, aes(x = long, y = lat))
map.occ <- map.occ + geom_tile(fill = 'cadetblue', colour = 'grey')
map.occ <- map.occ + geom_path(data = data_merged, 
                               mapping = aes(x = long, y = lat, 
                                             group = group, fill = NULL))
map.occ <- map.occ + coord_cartesian(xlim = c(-180, -20), ylim = c(5, 85))

