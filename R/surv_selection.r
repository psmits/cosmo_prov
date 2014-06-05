library(survival)
library(MuMIn)
library(stringr)
library(plyr)

source('surv_parametric.r')

var.imp <- function(models) {
  preds <- lapply(models, function(x) as.character(x$call$formula[-2])[-1])
  preds <- lapply(preds, function(x) unlist(str_split(x, ' \\+ ')))
  uni.pred <- unique(unlist(preds))

  wts <- Weights(laply(models, AICc))

  rel <- c()
  for(ii in seq(length(uni.pred))) {
    wh <- lapply(preds, function(x) any(x %in% uni.pred[ii]))
    rel[ii] <- sum(wts[unlist(wh)])
  }
  
  incepts <- grepl('[0-9]', uni.pred)

  rel <- rel[!incepts]
  rel <- cbind(data.frame(imp = rel), pred = uni.pred[!incepts])
  rel
}

na.imp <- var.imp(na.mod)
er.imp <- var.imp(er.mod)
