library(plyr)

source('../R/pa_mung.r')
source('../R/coverage.r')

source('http://bio.mq.edu.au/~jalroy/SQS-3-3.R')

stcov <- ddply(pa.mat, .(stage), summarize,
               uu = coverage(table(occurrence.genus_name)))

past <- split(pa.mat, pa.mat$stage)
subab <- list()
for (ii in seq(length(past))) {
  ab <- table(past[[ii]]$occurrence.genus_name)
  subab[[ii]] <- sqs(ab, q = stcov[ii, 2] - 0.06)[3]
}
names(subab) <- names(past)


# by dietary category
padiet <- split(pa.mat, pa.mat$comdiet)

dietcov <- llply(padiet, function(x) {
                 ddply(x, .(stage), summarize,
                       uu = coverage(table(occurrence.genus_name)))})

dietst <- lapply(padiet, function(x) {
                 split(x, x$stage)})

dietab <- list()
for (ii in seq(length(dietst))) {
  oo <- dietst[[ii]]
  uu <- dietcov[[ii]]
  dietsub <- list()
  for (jj in seq(length(oo))) {
    ab <- table(oo[[jj]]$occurrence.genus_name)
    dietsub[[jj]] <- sqs(ab, q = uu[jj, 2] - 0.1)[3]
  }
  names(dietsub) <- uu[, 1]
  dietab[[ii]] <- dietsub
}


# by locomotor category
paloco <- split(pa.mat, pa.mat$life_habit)

lococov <- llply(paloco, function(x) {
                 ddply(x, .(stage), summarize,
                       uu = coverage(table(occurrence.genus_name)))})

locost <- lapply(paloco, function(x) {
                 split(x, x$stage)})

locoab <- list()
for (ii in seq(length(locost))) {
  oo <- locost[[ii]]
  uu <- lococov[[ii]]
  locosub <- list()
  for (jj in seq(length(oo))) {
    ab <- table(oo[[jj]]$occurrence.genus_name)
    locosub[[jj]] <- sqs(ab, q = uu[jj, 2] - 0.1)[3]
  }
  names(locosub) <- uu[, 1]
  locoab[[ii]] <- locosub
}
