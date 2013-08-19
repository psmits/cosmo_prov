source('../src/multiclass_roc.r')

# various training and selection support functions
require(caret)
ctrl <- trainControl(#method = 'LOOCV',
                     method = 'repeatedcv',
                     #classProbs = TRUE,
                     number = 10,
                     repeats = 10)

nnetFuncs <- caretFuncs
nnetFuncs$fit <- function(x, y, first, last, ...) {
  nnet(x, y, ...)
}
nnetFuncs$summary <- multiClassSummary
nnet.ctrl <- rfeControl(functions = nnetFuncs,
                        method = 'repeatedcv',
                        repeats = 10,
                        number = 10,
                        verbose = FALSE,
                        returnResamp = 'final',
                        allowParallel = FALSE)

rfFuncs$summary <- multiClassSummary
rf.ctrl <- rfeControl(functions = rfFuncs,
                      method = 'repeatedcv',
                      repeats = 10,
                      number = 10,
                      verbose = FALSE, 
                      returnResamp = 'final')


make.form <- function(vari, resp) {
  form <- vector(mode = 'list', length = length(vari))
  for (ii in seq(length(vari))) {
    form[[ii]] <- as.formula(paste(paste(resp, ' ~ ', collapse = ''),
                                   paste(vari[seq(ii)],
                                         collapse = '+')))
  }
  form
}

make.form <- function(vari, resp) {
  form <- vector(mode = 'list', length = length(vari))
  for (ii in seq(length(vari))) {
    form[[ii]] <- as.formula(paste(paste(resp, ' ~ ', collapse = ''),
                                   paste(vari[seq(ii)],
                                         collapse = '+')))
  }
  form
}

data.maker <- function(gr, data, p = 0.75) {
  nd <- data[, colnames(data) %in% gr]
  out <- apply(nd, 2, createDataPartition,
               p = p, list = FALSE)
  out
}

multi.train <- function(form, data, seed = 1, ...) {
  # train across multiple formulas with the same data set
  #
  # Args
  #   form: list of formula objects
  #   data: data frame
  #   seed: random seed
  #   ...: arguments (matched exactly) for train.formula from caret
  #
  # Returns:
  #   list of model training results
  set.seed(seed)
  rr <- lapply(form, train, data = data, ...)
  rr
}


# plotting functions used in the paper

clean.resamp <- function(values) {
  s <- values
  s.m <- melt(s)
  nn <- data.frame(matrix(unlist(strsplit(as.character(s.m[, 2]), 
                                          split = '[~]', perl = TRUE)),
                          ncol = 2, byrow = TRUE))
  s.m <- cbind(s.m[, -2], nn)
  s.m <- split(s.m, f = s.m$X2)
  s.mod <- s.m$Accuracy$X1
  s.acc <- s.m$Accuracy$value
  s.kap <- s.m$Kappa$value
  s.dat <- data.frame(model = s.mod,
                      accuracy = s.acc,
                      kappa = s.kap)
  s.dat <- melt(s.dat)
  s.dat <- ddply(s.dat, .(model, variable), summarize,
                  mean = mean(value),
                  #se = sd(value)# / sqrt(length(value))
                  low = quantile(value, probs = 0.025),
                  high = quantile(value, probs = 0.975)
                  )
  s.dat
}

resamp.plot <- function(dat) {
  g.nre.s <- ggplot(dat, aes(x = model, y = mean))
  g.nre.s <- g.nre.s + geom_pointrange(mapping = aes(ymax = high,
                                                     ymin = low))
  g.nre.s <- g.nre.s + coord_flip()
  g.nre.s <- g.nre.s + facet_wrap(~ variable)
  g.nre.s
}

clean.diff <- function(diffs) {
  diffs <- melt(diffs)
  out <- ddply(diffs, .(X2, L1), summarize
               , mean = mean(value)
               #, se = sd(value)# / sqrt(length(value))
               #, low = t.test(value, conf.level = conf)$conf.int[1]
               #, high = t.test(value, conf.level = conf)$conf.int[2]
               , low = quantile(value, probs = 0.025)
               , high = quantile(value, probs = 0.975)
               )
  out
}

diffs.plot <- function(dat) {
  gdd <- ggplot(dat, aes(x = X2, y = mean))
  gdd <- gdd + geom_pointrange(mapping = aes(ymax = high,
                                             ymin = low))
  gdd <- gdd + geom_hline(yintercept = 0) + labs(x = 'Model')
  gdd <- gdd + coord_flip()
  gdd <- gdd + facet_wrap(~ L1)
  gdd
}

make.ptab <- function(val) {
  pp <- melt(lapply(val, #geo.rf.redi$stage$statistics$Accuracy, 
                      function(x) x$p.value))
  mods <- Reduce(rbind, strsplit(pp$L1, split = '[\\.]', perl = TRUE))
  pp <- cbind(as.data.frame(pp[, 1]), mods[, c(1, 3)])
  colnames(pp) <- c('p', 'mod1', 'mod2')
  ptable <- dcast(melt(pp, id.vars = c('mod1', 'mod2')), mod1 ~ mod2)
  ptab <- xtable(ptable)
}

grab.letters <- function(string, num) {
  str.vec <- strsplit(string, '')
  str.vec[[1]][seq(num)]
}

let2str <- function(x) {
  paste(x, collapse = '')
}

shorten <- function(x, n = 3) {
  let2str(grab.letters(x, n))
}
 

## insert copy statement here

## various functions to assist in the analysis accomplished here

summer <- function(cols, mat) {
  if (length(cols) != 0) {
    out <- rowSums(mat[, cols, drop = FALSE])
  } else {
    out <- rep(0, times = nrow(mat))
  }
  out
}


coverage <- function(x, dom = TRUE) {
  # calculate Good's U for count information
  #
  # Args:
  #  x: vector of abundances
  #  dom: include dominate occurence? (default TRUE)
  #
  # Returns:
  #  value between 0 and 1 corresponding to Good's U

  O <- ifelse(dom, sum(x), (sum(x) - max(x)))
  n <- sum(x == 1)

  u <- 1 - (n / O)

  u
}

eco.covmin <- function(abund, cats) {
  # calculate the min coverage for some category from abundance matrix
  #
  # Args:
  #  abund: an abundance matrix
  #  cats: a vector of categories
  #
  # Returns:
  #  list of minimum coverage values for each category
  #  names of list are category

  ulist <- by(abund, cats, function(x) apply(x, 2, coverage), 
              simplify = FALSE)
  ulist <- lapply(ulist,  function(x) ifelse(is.nan(x), NA, x))
  umin <- lapply(ulist, function(x) min(x[x > 0], na.rm = TRUE))
  umin
}

eco.bin <- function(abund, cats) {
  # sum the total occurences for some categories
  #
  # Args:
  #  abund: a vector of abundances
  #  cats: a vector of categories
  #
  # Returns:
  #  vector of category abundance with names

  suppressMessages(dats <- melt(data.frame(abund, cats))[, c(1, 3)])
  ## find all of one type
  ## sum that
  uni <- unique(cats)
  cat.abund <- array(dim = length(uni))
  names(cat.abund) <- uni
  for (ii in seq(length(uni))) {
    cat.abund[ii] <- sum(dats[dats[, 1] == uni[ii], 2])
  }

  cat.abund
}

bin.sub <- function(abund, cats, u, skip = TRUE) {
  # THIS FUNCTION DOES NOT WORK
  # THERE IS A PROBLEM WITH SQS
  # need to write in a skipper for things that can never reach that coverage
  # cut up some factor by a factor then subsample them at a coverage
  #
  # Args:
  #  abund: vector of abundance information
  #  cats: vector of categories. length == length(abund)
  #  u: minimum coverage for each cat. length == unique(cats)
  #  skip: skip if cats u is below global u. default = TRUE
  #
  # Returns:
  #  mean subsampled richness for each category

  cut.up <- split(abund, cats)
  sub.abund <- vector(mode = 'list', length = length(cut.up))
  names(sub.abund) <- names(cut.up)

  for (ii in seq(length(u))) {
    ## need to place this in a list for each category
    ## currently just only produces the value for the last one
    if (sum(cut.up[[ii]]) == 0) {
      sub.abund[[ii]] <- 0
    } else if (coverage(cut.up[[ii]]) < u[ii]) {
      sub.abund[[ii]] <- NA
    } else {
      sub.abund[[ii]] <- sqs(cut.up[[ii]], q = u[ii])['subsampled richness']
    }
  }

  sub.abund
}

