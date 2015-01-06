library(rstan)
library(arm)
library(parallel)
library(xtable)
library(stringr)
library(plyr)
library(igraph)
library(ape)


RNGkind(kind = "L'Ecuyer-CMRG")
seed <- 420
nsim <- 1000

set.seed(seed)
source('../R/surv_setup.r')
load('../data/taxonomy_tree.rdata')

# compile the models
poisson.mod <- stan(file = '../stan/degree_model.stan') # 
pois.phy.mod <- stan(file = '../stan/degree_phy_model.stan') # 
pois.spt.mod <- stan(file = '../stan/degree_spt_model.stan') # 
pois.ful.mod <- stan(file = '../stan/degree_full_model.stan') # 

negbin.mod <- stan(file = '../stan/deg_mod_over.stan')
negb.phy.mod <- stan(file = '../stan/deg_phy_over.stan')
negb.phy.mod <- stan(file = '../stan/deg_spt_over.stan')
negb.ful.mod <- stan(file = '../stan/deg_full_over.stan')


# prep the data
species.graph <- function(bipartite) {
  bip <- bipartite.projection(bipartite)
  if(any(grepl('[0-9]', V(bip[[1]])$name))) {
    taxa <- bip[[2]]
    st <- bip[[1]] 
  } else {
    taxa <- bip[[1]]
    st <- bip[[2]] 
  }
  out <- taxa
  out
}

spg <- llply(taxawin, function(x) {
             o <- species.graph(x)
             E(o)$weight <- 1
             o})
name <- llply(spg, function(x) V(x)$name)
inin <- llply(name, function(x) x %in% na.ecol$taxa)
spg <- Map(function(x, y) delete.vertices(x, which(!y)), spg, inin)

adj <- llply(spg, function(x) get.adjacency(x, sparse = FALSE))
deg <- llply(spg, degree)
size <- llply(spg, vcount)
name <- llply(spg, function(x) V(x)$name)

# have to make these line up correctly
ecol <- llply(name, function(x) {
              tt <- na.ecol[na.ecol$taxa %in% x, ]
              tt[match(x, tt[, 1]), ]
              tt})
ecol <- llply(ecol, function(x) {
              x$mass <- rescale(log(x$mass))
              x})

diet <- llply(ecol, function(x) model.matrix( ~ x$diet - 1)[, -1])
move <- llply(ecol, function(x) model.matrix( ~ x$move - 1)[, -1])


# prepare the phylo vcv
tax.nam <- str_replace(na.ecol[, 1], ' ', '_')
to.drop <- na.scale$tip.label[!(na.scale$tip.label %in% tax.nam)]
na.tree <- drop.tip(na.scale, to.drop)

# for each window
vcv.s <- list()
for(ii in seq(length(name))) {
  nono <- na.tree$tip.label[!(na.tree$tip.label %in% 
                              str_replace(name[[ii]], ' ', '_'))]
  temp <- drop.tip(na.tree, nono)
  temp$edge.length <- temp$edge.length / max(diag(vcv(temp)))
  temp.vcv <- vcv(temp)
  
  # have to have this line up correctly
  o <- match(str_replace(name[[ii]], ' ', '_'), colnames(temp.vcv))
  temp.vcv <- temp.vcv[o, o]

  vcv.s[[ii]] <- temp.vcv  
}


# fit the models
# basic poisson model
degfit <- list()
for(ii in seq(length(adj))) {

  data <- list(N = size[[ii]],
               D = ncol(diet[[ii]]),
               M = ncol(move[[ii]]),
               degree = deg[[ii]],
               mass = ecol[[ii]]$mass,
               diet = diet[[ii]],
               move = move[[ii]],
               vcv = vcv.s[[ii]],
               adj = adj[[ii]])

  degmod <- mclapply(1:4, mc.cores = detectCores(),
                     function(x) stan(fit = poisson.mod, 
                                      seed = seed,
                                      data = data, 
                                      chains = 1, chain_id = x,
                                      refresh = -1))
  degfit[[ii]] <- sflist2stanfit(degmod)
}





# probably need to switch to a negative binomial, duh :P
overfit <- list()
for(ii in seq(length(adj))) {
  data <- list(N = size[[ii]],
               D = ncol(diet[[ii]]),
               M = ncol(move[[ii]]),
               degree = deg[[ii]],
               mass = ecol[[ii]]$mass,
               diet = diet[[ii]],
               move = move[[ii]],
               vcv = vcv.s[[ii]],
               adj = adj[[ii]])

  degmod <- mclapply(1:4, mc.cores = detectCores(),
                     function(x) stan(fit = negbin.mod, 
                                      seed = seed,
                                      data = data, 
                                      chains = 1, chain_id = x,
                                      refresh = -1))
  overfit[[ii]] <- sflist2stanfit(degmod)
}
