library(ggplot2)
library(reshape2)
library(scales)
library(hexbin)
library(stringr)
library(grid)
library(survival)
library(GGally)
library(moments)
library(plyr)

nsim <- 100

pairwise.diffs <- function(x) {
  # create column combination pairs
  prs <- cbind(rep(1:ncol(x), each = ncol(x)), 1:ncol(x))
  col.diffs <- prs[prs[, 1] < prs[, 2], , drop = FALSE]
  # do pairwise differences 
  result <- x[, col.diffs[, 1]] - x[, col.diffs[, 2], drop = FALSE]
  # set colnames
  if(is.null(colnames(x)))
    colnames(x) <- 1:ncol(x)
  
  colnames(result) <- paste(colnames(x)[col.diffs[, 1]], ".vs.", 
                            colnames(x)[col.diffs[, 2]], sep = "")
  result
}



theme_set(theme_bw())
cbp <- c('#E69F00', '#56B4E9', '#009E73', '#F0E442', 
         '#0072B2', '#D55E00', '#CC79A7')
theme_update(axis.text = element_text(size = 20),
             axis.title = element_text(size = 30),
             legend.text = element_text(size = 25),
             legend.title = element_text(size = 26),
             legend.key.size = unit(2, 'cm'),
             strip.text = element_text(size = 20))

# discrepency
# on point
pnt.est <- llply(negbinom.cout, function(y) 
                 aaply(y, .margins = 2, .fun = median))
pnt.est <- cbind(melt(Reduce(rbind, pnt.est)), label = rep(1:test, 3))[, -1]

mgt.est <- llply(negbinom.cout, function(y) 
                 aaply(y, .margins = 2, .fun = function(x) 
                       quantile(x, c(0.2, 0.8))))
mgt.est <- cbind(Reduce(rbind, mgt.est), label = rep(1:test, each = 3))
mgt.est <- cbind(data.frame(mgt.est), quant = rownames(mgt.est))
names(mgt.est) <- c('bot', 'top', 'label', 'quant')
mgt.est <- Reduce(rbind, llply(split(mgt.est, mgt.est$quant), 
                               function(x) x[order(x$label), ]))
pnt.est <- cbind(pnt.est, mgt.est[, 1:2])
names(pnt.est) <- c('quant', 'value', 'label', 'bot', 'top')
pnt.est <- pnt.est[!(abs(pnt.est$top - pnt.est$bot) > 100), ]

bad.one <- llply(split(pnt.est, pnt.est$quant), 
                 function(x) which(!(1:10 %in% x$label)))
bad.one <- melt(bad.one)
names(bad.one) <- c('label', 'quant')
bad.one$value <- rep(0, nrow(bad.one))
bad.one$label <- as.numeric(as.character(bad.one$label))

# forward
pnt.for <- llply(negbinom.cfor, function(y)
                 aaply(y, .margins = 2, .fun = median))
pnt.for <- cbind(melt(Reduce(rbind, pnt.for)), 
                 label = rep(1:(test - 1), 3))[, -1]

mgt.for <- llply(negbinom.cfor, function(y) 
                 aaply(y, .margins = 2, .fun = function(x) 
                       quantile(x, c(0.2, 0.8))))
mgt.for <- cbind(Reduce(rbind, mgt.for), 
                 label = rep(1:(test - 1), each = 3))
mgt.for <- cbind(data.frame(mgt.for), quant = rownames(mgt.for))
names(mgt.for) <- c('bot', 'top', 'label', 'quant')
mgt.for <- Reduce(rbind, llply(split(mgt.for, mgt.for$quant), 
                               function(x) x[order(x$label), ]))
pnt.for <- cbind(pnt.for, mgt.for[, 1:2])
names(pnt.for) <- c('quant', 'value', 'label', 'bot', 'top')
pnt.for <- pnt.for[!(abs(pnt.for$top - pnt.for$bot) > 100),]

bad.two <- llply(split(pnt.for, pnt.for$quant), 
                 function(x) which(!(1:10 %in% x$label)))
bad.two <- melt(bad.two)
names(bad.two) <- c('label', 'quant')
bad.two$value <- rep(0, nrow(bad.two))
bad.two$label <- as.numeric(as.character(bad.two$label))

pnt.est <- rbind(cbind(pnt.est, mod = rep('est', nrow(pnt.est))),
                 cbind(pnt.for, mod = rep('for', nrow(pnt.for))))

bad.one <- rbind(cbind(bad.one, mod = rep('est', nrow(bad.one))),
                 cbind(bad.two, mod = rep('for', nrow(bad.two))))

disc <- ggplot(pnt.est, aes(x = label, y = value, 
                            ymin = bot, ymax = top,
                            colour = mod))
disc <- disc + geom_hline(aes(yintercept = 0), colour = 'grey', size = 2)
disc <- disc + geom_pointrange(position = position_jitter(height = 0, 
                                                          width = 0.1))
disc <- disc + geom_point(data = bad.one, 
                          mapping = aes(x = label, y = value, 
                                        ymin = value, ymax = value),
                          size = 5, shape = 8, 
                          position = position_jitter(height = 0, 
                                                     width = 0.1))
disc <- disc + facet_grid(quant ~ .)
disc <- disc + scale_colour_manual(values = cbp)
