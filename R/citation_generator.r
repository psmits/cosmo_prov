library(knitcitations)

libs <- c('RCore' = citation(),
          'igraph' = citation('igraph'),
          'treebase' = citation('treebase'),
          'rstan' = citation('rstan'),
          'taxize' = citation('taxize'),
          'ape' = citation('ape'),
          'phytools' = citation('phytools'))



write.bibtex(libs, file = '../doc/packages.bib')
