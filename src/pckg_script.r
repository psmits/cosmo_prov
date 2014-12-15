options(repos = 'http://cran.r-project.org')

Sys.setenv(MAKEFLAGS = "-j4")
source("http://mc-stan.org/rstan/install.R") 
install_rstan()

packs <- list('arm', 'ape', 'parallel', 'stringr', 'igraph', 'maps', 
              'plyr', 'survival', 'reshape2', 'mapproj', 'phytools',
              'paleotree', 'taxize', 'devtools')
lapply(packs, function(x)  
       install.packages(x, dependencies = TRUE, 
                        repos = 'http://cran.r-project.org'))

q()
