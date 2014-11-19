library(treebase)
library(parallel)

RNGkind(kind = "L'Ecuyer-CMRG")
seed <- 420

source('../R/na_mung.r')


# get as many trees as possible
tax <- unique(dat$occurrence.genus_name)

tax.out <- mclapply(tax, function(x) 
                    try(search_treebase(x, returns = 'tree', 
                                        by = 'taxon', 
                                        verbose = FALSE, 
                                        only_metadata = TRUE)), 
                    mc.cores = detectCores())

out <- tax.out[!str_detect(laply(tax.out, class), 'try')]
flat <- unlist(out, recursive = FALSE)
sp <- laply(flat, function(x) x$kind)
flat <- flat[which(sp == 'Species Tree')]
flat <- laply(flat, function(x) x$Tr.id)
flat <- unique(flat)

mam.trees <- mclapply(flat, function(x) 
                    try(search_treebase(x, returns = 'tree', 
                                        by = 'id.tree', 
                                        verbose = FALSE)), 
                    mc.cores = detectCores())


mam.flat <- unlist(mam.trees, recursive = FALSE)
class(mam.flat) <- 'multiPhylo'
save(mam.flat, file = '../data/many_trees.rdata')
