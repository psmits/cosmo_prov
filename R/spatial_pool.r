library(rstan)
library(arm)
library(parallel)
library(stringr)
library(plyr)
library(igraph)
library(ape)

RNGkind(kind = "L'Ecuyer-CMRG")
seed <- 420
nsim <- 1000

set.seed(seed)
source('../R/na_mung.r')
load('../data/scaled_super.rdata')

adj <- matrix(0, nrow = ncell(sp.ras), ncol = ncell(sp.ras))
adj.list <- adjacent(sp.ras, seq(ncell(sp.ras)), directions = 8)
adj.list <- unique(adj.list)
by.cell <- split(adj.list, adj.list[, 1])
keeps <- Map(function(x, y) !(x %in% as.numeric(y)), 
             x = by.cell, y = names(by.cell))
neigh <- laply(Map(function(x, y) x[y], x = by.cell, y = keeps), length)
adj[adj.list] <- 1

splits <- split(dat, dat$bins)
diversities <- llply(splits, function(x) table(x$gid))


data <- list()
for(ii in seq(length(diversities))) {
  div <- array(data = 0, dim = nrow(adj))
  div[as.numeric(names(diversities[[ii]]))] <- diversities[[ii]]

  data[[ii]] <- list(N = length(div), diversity = div, 
                    adj = adj, neighbors = neigh)

  nam <- paste0('../data/data_dump/spatial_info_', ii, '.data.R')
  with(data[[ii]], {stan_rdump(list = c('N', 'diversity', 'adj', 'neighbors'),
                               file = nam)})
}

#par(mfrow = c(6, 5), mar = c(4, 4, 2, 2))
#for(ii in seq(30)) {
#  values(sp.ras) <- as.vector(data[[ii]]$div)
#  plot(sp.ras, main = ii)
#}

#fit.mod <- stan(file = '../stan/spatial_model.stan')
#fit.list <- mclapply(1:4, mc.cores = detectCores(),
#                     function(x) stan(fit = fit.mod,
#                                      data = data[[1]],
#                                      seed = 420,
#                                      chains = 1, chain_id = x,
#                                      refresh = -1))
#fit.attempt <- sflist2stanfit(fit.list)
