library(ape)
library(plyr)
library(reshape2)
library(ggplot2)
library(stringr)

source('../R/surv_setup.r')
load('../data/scaled_super.rdata')
load('../data/update_taxonomy.rdata')

# north america
no.fam <- which(na.tax$family_name == '')
no.ord <- which(na.tax$order_name == '')
rms <- unique(c(no.fam, no.ord))
clean.na <- na.tax[-rms, ]
uni.tax <- unique(clean.na[, 1:4])

new.tax <- replace.taxonomy(dat, uni.tax[, 1:3])
new.tax$order_name[new.tax$order_name == ''] <- 'Unk'
new.tax$family_name[new.tax$family_name == ''] <- 'Unk'
for(ii in seq(nrow(new.tax))) {
  if(new.tax$order_name[ii] == 'Artiodactyla') {
    new.tax$order_name[ii] <- 'Cetartiodactyla'
  }
  if(new.tax$family_name[ii] == 'Unk') {
    new.tax$family_name[ii] <- paste0(new.tax$order_name[ii], 
                                      new.tax$family_name[ii])
  }
}

my.taxonomy <- new.tax[, c('order_name', 'family_name', 
                           'occurrence.genus_name', 'name.bi')]
my.taxonomy <- unique(my.taxonomy)
mixx <- str_replace(my.taxonomy$name.bi, ' ', '_')

re.tax <- my.taxonomy[match(spt$tip.label, mixx), ]
re.tax <- data.frame(re.tax)

names(re.tax)[4] <- 'taxa'
row.names(re.tax) <- NULL
#length(unique(re.tax$order_name))


plot.phylo(spt, show.tip.label = FALSE)
#tiplabels(pch = 21, bg = 'blue', col = 'blue', cex = .5)
