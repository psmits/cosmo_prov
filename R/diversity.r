  library(plyr)

  source('../R/na_mung.r')
  source('../R/coverage.r')

  source('http://bio.mq.edu.au/~jalroy/SQS-3-3.R')

  stcov <- ddply(pa.mat, .(stage), summarize,
                 uu = coverage(table(occurrence.genus_name)))

  past <- split(pa.mat, pa.mat$stage)
  subab <- list()
  for (ii in seq(length(past))) {
    ab <- table(past[[ii]]$occurrence.genus_name)
    subab[[ii]] <- sqs(ab, q = min(stcov[, 2]) - 0.06)[3]
  }
  names(subab) <- names(past)

  # do at 2My bins
  bins <- seq(from = 0, to = 66, by = 2)
  bins <- cbind(top = bins[-1], bot = bins[-length(bins)])
  # assign every occurence to a bin
  pa.mat$bins <- rep(NA, nrow(pa.mat))
  for (ii in seq(nrow(bins))) {
    out <- which(pa.mat$ma_mid < bins[ii, 1] & pa.mat$ma_mid >= bins[ii, 2])
    pa.mat$bins[out] <- bins[ii, 1]
  }
  bincov <- ddply(pa.mat, .(bins), summarize,
                  uu = coverage(table(occurrence.genus_name)))

  pabin <- split(pa.mat, pa.mat$bins)
  binsub <- list()
  for (ii in seq(length(pabin))) {
    ab <- table(pabin[[ii]]$occurrence.genus_name)
    binsub[[ii]] <- sqs(ab, q = min(bincov[, 2]) - 0.06)[3]
  }


  # by dietary category
  padiet <- split(pa.mat, pa.mat$comdiet)

  dietcov <- llply(padiet, function(x) {
                   ddply(x, .(stage), summarize,
                         uu = coverage(table(occurrence.genus_name)))})

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
      dietsub[[jj]] <- sqs(ab, q = min(uu[, 2]) - 0.1)[3]
    }
    names(dietsub) <- uu[, 1]
    dietab[[ii]] <- dietsub
  }
  names(dietab) <- names(dietst)

  # repeate at bin level
  dtbncov <- llply(padiet, function(x) {
                   ddply(x, .(bins), summarize,
                       uu = coverage(table(occurrence.genus_name)))})

dtbin <- lapply(padiet, function(x) {
                split(x, x$bins)})
dtbinab <- list()
for (ii in seq(length(dtbin))) {
  oo <- dtbin[[ii]]
  uu <- dtbncov[[ii]]
  dtbnsb <- list()
  for (jj in seq(length(oo))) {
    ab <- table(oo[[jj]]$occurrence.genus_name)
    dtbnsb[[jj]] <- sqs(ab, q = min(uu[, 2]) - 0.1)[3]
  }
  names(dtbnsb) <- uu[, 1]
  dtbinab[[ii]] <- dtbnsb
}
names(dtbinab) <- names(dtbin)


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
    locosub[[jj]] <- sqs(ab, q = min(uu[, 2]) - 0.1)[3]
  }
  names(locosub) <- uu[, 1]
  locoab[[ii]] <- locosub
}
names(locoab) <- names(locost)
