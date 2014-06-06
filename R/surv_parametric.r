library(survival)
library(MuMIn)
source('../R/model_sel.r')

source('../R/na_surv.r')
source('../R/eur_surv.r')

na.ecol$mass <- log(as.numeric(as.character(na.ecol$mass)))
#er.ecol$mass <- log(as.numeric(as.character(er.ecol$mass)))
na.genecol$mass <- log(as.numeric(as.character(na.genecol$mass)))
#er.genecol$mass <- log(as.numeric(as.character(er.genecol$mass)))

nsim <- 1000

# just north america
# species
vars <- names(na.ecol)[-1]
mods <- create.model(vars = vars)

surv.wei <- survreg(na.surv ~ 1, data = na.ecol, dist = 'weibull')
na.wei <- fit.models(surv.wei, mods, data = na.ecol, dist = 'weibull')

pna.wei <- list()
for(ii in seq(nsim)) {
  pna.wei[[ii]] <- var.imp(model.base(na.surv, mods = mods, 
                                      datas = na.ecol, 
                                      distribution = 'weibull'))
}

surv.exp <- survreg(na.surv ~ 1, data = na.ecol, dist = 'exponential')
na.exp <- fit.models(surv.exp, mods, data = na.ecol, dist = 'exponential')

pna.exp <- list()
for(ii in seq(nsim)) {
  pna.exp[[ii]] <- var.imp(model.base(na.surv, mods = mods, 
                                      datas = na.ecol, 
                                      distribution = 'weibull'))
}

na.mod <- c(na.wei, na.exp)
pna.mod <- c(pna.wei, pna.exp)

# genera
vars <- names(na.genecol)[-c(1)]
mods <- create.model(vars = vars)

surv.wei <- survreg(nagen.surv ~ 1, data = na.genecol, dist = 'weibull')
nagen.wei <- fit.models(surv.wei, mods, data = na.genecol, dist = 'weibull')

pnagen.wei <- list()
for(ii in seq(nsim)) {
  pnagen.wei[[ii]] <- var.imp(model.base(nagen.surv, mods = mods, 
                                         datas = na.genecol, 
                                         distribution = 'weibull'))
}

surv.exp <- survreg(nagen.surv ~ 1, data = na.genecol, dist = 'exponential')
nagen.exp <- fit.models(surv.wei, mods, data = na.genecol, dist = 'exponential')

pnagen.exp <- list()
for(ii in seq(nsim)) {
  pnagen.exp[[ii]] <- var.imp(model.base(nagen.surv, mods = mods, 
                                         datas = na.genecol, 
                                         distribution = 'exponential'))
}

nagen.mod <- c(nagen.wei, nagen.exp)
pnagen.mod <- c(pnagen.wei, pnagen.exp)


# just europe
# species
vars <- names(er.ecol)[-1]
mods <- create.model(vars = vars)

surv.wei <- survreg(er.surv ~ 1, data = er.ecol, dist = 'weibull')
er.wei <- fit.models(surv.wei, mods, data = er.ecol, dist = 'weibull')

per.wei <- list()
for(ii in seq(nsim)) {
  per.wei[[ii]] <- var.imp(model.base(er.surv, mods = mods, 
                                      datas = er.ecol, 
                                      distribution = 'weibull'))
}

surv.exp <- survreg(er.surv ~ 1, data = er.ecol, dist = 'exponential')
er.exp <- fit.models(surv.exp, mods, data = er.ecol, dist = 'exponential')

per.exp <- list()
for(ii in seq(nsim)) {
  per.exp[[ii]] <- var.imp(model.base(er.surv, mods = mods, 
                                      datas = er.ecol, 
                                      distribution = 'exponential'))
}

er.mod <- c(er.wei, er.exp)
per.mod <- c(per.wei, per.exp)


# genera
vars <- names(er.genecol)[-1]
mods <- create.model(vars = vars)

surv.wei <- survreg(ergen.surv ~ 1, data = er.genecol, dist = 'weibull')
ergen.wei <- fit.models(surv.wei, mods, data = er.genecol, dist = 'weibull')

pergen.wei <- list()
for(ii in seq(nsim)) {
  pergen.wei[[ii]] <- var.imp(model.base(ergen.surv, mods = mods, 
                                         datas = er.genecol, 
                                         distribution = 'weibull'))
}

surv.exp <- survreg(ergen.surv ~ 1, data = er.genecol, dist = 'exponential')
ergen.exp <- fit.models(surv.exp, mods, data = er.genecol, dist = 'exponential')

pergen.exp <- list()
for(ii in seq(nsim)) {
  pergen.exp[[ii]] <- var.imp(model.base(ergen.surv, mods = mods, 
                                         datas = er.genecol, 
                                         distribution = 'exponential'))
}

ergen.mod <- c(ergen.wei, ergen.exp)
pergen.mod <- c(pergen.wei, pergen.exp)
