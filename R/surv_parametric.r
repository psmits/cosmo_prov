library(survival)
library(MuMIn)

source('../R/na_surv.r')
source('../R/eur_surv.r')

#na.ecol$mass <- log(as.numeric(as.character(na.ecol$mass)))
#er.ecol$mass <- log(as.numeric(as.character(er.ecol$mass)))


# just north america
# species
vars <- names(na.ecol)[-1]
mod <- list()
nll <- na.surv ~ 1
for(ii in seq(length(vars))) {
  cc <- combn(vars, ii, simplify = FALSE)
  mod[[ii]] <- cc #lapply(cc, function(x) paste(x, collapse = ' + '))
}
mod <- unlist(mod, recursive = FALSE)

surv.wei <- survreg(na.surv ~ 1, data = na.ecol, dist = 'weibull')
ups <- list()
ups[[1]] <- surv.wei
for(ii in seq(from = 2, to = (length(mod) + 1))) {
  ups[[ii]] <- update(surv.wei, paste('. ~ . +', 
                                      paste(mod[[ii - 1]], collapse = ' + ')))
}
na.wei <- ups

surv.exp <- survreg(na.surv ~ 1, data = na.ecol, dist = 'exponential')
ups <- list()
ups[[1]] <- surv.exp
for(ii in seq(from = 2, to = (length(mod) + 1))) {
  ups[[ii]] <- update(surv.exp, paste('. ~ . +', 
                                      paste(mod[[ii - 1]], collapse = ' + ')))
}
na.exp <- ups

na.mod <- c(na.wei, na.exp)

# just europe
# species
vars <- names(er.ecol)[-1]
mod <- list()
for(ii in seq(length(vars))) {
  cc <- combn(vars, ii, simplify = FALSE)
  mod[[ii]] <- cc #lapply(cc, function(x) paste(x, collapse = ' + '))
}
mod <- unlist(mod, recursive = FALSE)

surv.wei <- survreg(er.surv ~ 1, data = er.ecol, dist = 'weibull')
ups <- list()
ups[[1]] <- surv.wei
for(ii in seq(from = 2, to = (length(mod) + 1))) {
  ups[[ii]] <- update(surv.wei, paste('. ~ . +', 
                                      paste(mod[[ii - 1]], collapse = ' + ')))
}
er.wei <- ups

surv.exp <- survreg(er.surv ~ 1, data = er.ecol, dist = 'exponential')
ups <- list()
ups[[1]] <- surv.exp
for(ii in seq(from = 2, to = (length(mod) + 1))) {
  ups[[ii]] <- update(surv.exp, paste('. ~ . +', 
                                      paste(mod[[ii - 1]], collapse = ' + ')))
}
er.exp <- ups

er.mod <- c(er.wei, er.exp)
