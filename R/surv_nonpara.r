library(survival)
library(MuMIn)

source('../R/na_surv.r')
source('../R/eur_surv.r')

na.ecol$mass <- log(as.numeric(as.character(na.ecol$mass)))
#er.ecol$mass <- log(as.numeric(as.character(er.ecol$mass)))

# north america
# species
vars <- names(na.ecol)[-1]
mod <- list()
nll <- na.surv ~ 1
for(ii in seq(length(vars))) {
  cc <- combn(vars, ii, simplify = FALSE)
  mod[[ii]] <- cc #lapply(cc, function(x) paste(x, collapse = ' + '))
}
mod <- unlist(mod, recursive = FALSE)

nakm <- survfit(formula = na.surv ~ 1, data = na.ecol)
ups <- list()
ups[[1]] <- nakm
for(ii in seq(from = 2, to = (length(mod) + 1))) {
  ups[[ii]] <- update(nakm, paste('. ~ . +', 
                                  paste(mod[[ii - 1]], collapse = ' + ')))
}
na.species <- ups

# genera
vars <- names(na.genecol)[-1]
mod <- list()
nll <- nagen.surv ~ 1
for(ii in seq(length(vars))) {
  cc <- combn(vars, ii, simplify = FALSE)
  mod[[ii]] <- cc #lapply(cc, function(x) paste(x, collapse = ' + '))
}
mod <- unlist(mod, recursive = FALSE)

nakm <- survfit(formula = nagen.surv ~ 1, data = na.genecol)
ups <- list()
ups[[1]] <- nakm
for(ii in seq(from = 2, to = (length(mod) + 1))) {
  ups[[ii]] <- update(nakm, paste('. ~ . +', 
                                  paste(mod[[ii - 1]], collapse = ' + ')))
}
na.genera <- ups


# europe
# species
vars <- names(er.ecol)[-1]
mod <- list()
nll <- er.surv ~ 1
for(ii in seq(length(vars))) {
  cc <- combn(vars, ii, simplify = FALSE)
  mod[[ii]] <- cc #lapply(cc, function(x) paste(x, collapse = ' + '))
}
mod <- unlist(mod, recursive = FALSE)

erkm <- survfit(formula = er.surv ~ 1, data = er.ecol)
ups <- list()
ups[[1]] <- erkm
for(ii in seq(from = 2, to = (length(mod) + 1))) {
  ups[[ii]] <- update(erkm, paste('. ~ . +', 
                                  paste(mod[[ii - 1]], collapse = ' + ')))
}
er.species <- ups

# genera
vars <- names(er.genecol)[-1]
mod <- list()
nll <- ergen.surv ~ 1
for(ii in seq(length(vars))) {
  cc <- combn(vars, ii, simplify = FALSE)
  mod[[ii]] <- cc #lapply(cc, function(x) paste(x, collapse = ' + '))
}
mod <- unlist(mod, recursive = FALSE)

erkm <- survfit(formula = ergen.surv ~ 1, data = er.genecol)
ups <- list()
ups[[1]] <- erkm
for(ii in seq(from = 2, to = (length(mod) + 1))) {
  ups[[ii]] <- update(erkm, paste('. ~ . +', 
                                  paste(mod[[ii - 1]], collapse = ' + ')))
}
er.genera <- ups
