options(repos = 'http://cran.r-project.org')

Sys.setenv(MAKEFLAGS = "-j4")
source("http://mc-stan.org/rstan/install.R") 
install_rstan()

packs <- list('arm', 'ape', 'parallel', 'stringr', 'igraph', 'maps', 
              'plyr', 'survival', 'reshape2', 'mapproj')
install.packages(packs, dependencies = TRUE, 
                 repos = 'http://cran.r-project.org')
