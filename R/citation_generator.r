library(knitcitations)

libs <- list('RCore' = citation(),
             'paleoTS' = citation('paleoTS'),
             'igraph' = citation('igraph'),
             'boot' = citation('boot')[1])

write.bibtex(libs, file = '../doc/packages.bib')
