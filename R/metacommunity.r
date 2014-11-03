library(igraph)
library(plyr)
library(parallel)
library(rstan)
library(boot)

source('../R/cosmo_prov.r')

RNGkind(kind = "L'Ecuyer-CMRG")
seed <- 420

split.graph <- function(graph) {
  bip <- bipartite.projection(graph)

  if(any(grepl('[0-9]', V(bip[[1]])$name))) {
    taxa <- V(bip[[2]])$name
    st <- V(bip[[1]])$name 
  } else {
    taxa <- V(bip[[1]])$name
    st <- V(bip[[2]])$name 
  }

  out <- list()
  out$taxa <- taxa
  out$locality <- st
  out
}

node.taxa <- function(graph, data) {
  ex <- split.graph(graph)
  ntaxa <- length(ex$taxa)
  nloc <- length(ex$locality)

  in.bin <- data[data$name.bi %in% ex$taxa,]
  in.bin <- in.bin[!duplicated(in.bin$name.bi), ]
  in.bin <- in.bin[, c('name.bi', 'comdiet', 'comlife')]
  dd <- model.matrix( ~ comdiet - 1, in.bin)
  colnames(dd) <- c('carnivore', 'herbivore', 'insectivore', 'omnivore')
  ll <- model.matrix( ~ comlife - 1, in.bin)
  colnames(ll) <- c('arboreal', 'ground', 'scansorial')
  inter <- interaction(in.bin$comdiet, in.bin$comlife)
  inter <- model.matrix( ~ inter - 1)
  colnames(inter) <- c('carnivore.arboreal', 'herbivore.arboreal', 
                       'insectivore.arboreal', 'omnivore.arboreal',
                       'carnivore.ground', 'herbivore.ground',
                       'insectivore.ground', 'omnivore.ground',
                       'carnivore.scansorial', 'herbivore.scansorial',
                       'insectivore.scansorial', 'omnivore.scansorial')

  node.properties <- cbind(dd, ll, inter)
  node.properties
}

dat$comdiet <- factor(dat$comdiet)
dat$comlife <- factor(dat$comlife)
node.matrix <- lapply(taxawin, node.taxa, data = dat)
links <- laply(taxawin, ecount)
ntaxa <- laply(node.matrix, nrow)
nloc <- laply(taxawin, function(x) length(split.graph(x)$locality))
code.length <- unlist(win.bg$code)

preds <- laply(node.matrix, function(x) apply(x, 2, sum))
cats <- ncol(preds)

data <- list(T = length(code.length),
             N = ntaxa,
             L = nloc,
             C = code.length,
             E = links,
             np = cats,
             preds = preds)

model <- stan(file = '../stan/graph_time.stan')

modlist <- mclapply(1:4, mc.cores = detectCores(),
                    function(x) stan(fit = model, 
                                     seed = seed,
                                     data = data,
                                     chains = 1, chain_id = x,
                                     refresh = -1))

initmod <- sflist2stanfit(modlist)
