library(survival)
library(MuMIn)

source('../R/na_surv.r')
source('../R/eur_surv.r')

# just north america
# species
na.wei <- survreg(formula = na.surv ~ 1, dist = 'weibull')
nad.wei <- survreg(formula = na.surv ~ diet, data = na.ecol, dist = 'weibull')
nal.wei <- survreg(formula = na.surv ~ move, data = na.ecol, dist = 'weibull')
nadl.wei <- survreg(formula = na.surv ~ diet + move,
                    data = na.ecol, dist = 'weibull')
na.exp <- survreg(formula = na.surv ~ 1, dist = 'exponential')
nad.exp <- survreg(formula = na.surv ~ diet, data = na.ecol, dist = 'exponential')
nal.exp <- survreg(formula = na.surv ~ move, data = na.ecol, dist = 'exponential')
nadl.exp <- survreg(formula = na.surv ~ diet + move,
                    data = na.ecol, dist = 'exponential')

na.mod <- list(na.wei, nad.wei, nal.wei, nadl.wei,
               na.exp, nad.exp, nal.exp, nadl.exp)

# genus
nag.wei <- survreg(formula = nagen.surv ~ 1, dist = 'weibull')
nagd.wei <- survreg(formula = nagen.surv ~ diet, data = na.genecol, dist = 'weibull')
nagl.wei <- survreg(formula = nagen.surv ~ move, data = na.genecol, dist = 'weibull')
nagdl.wei <- survreg(formula = nagen.surv ~ diet + move,
                     data = na.genecol, dist = 'weibull')
nag.exp <- survreg(formula = nagen.surv ~ 1, dist = 'exponential')
nagd.exp <- survreg(formula = nagen.surv ~ diet, data = na.genecol, dist = 'exponential')
nagl.exp <- survreg(formula = nagen.surv ~ move, data = na.genecol, dist = 'exponential')
nagdl.exp <- survreg(formula = nagen.surv ~ diet + move,
                     data = na.genecol, dist = 'exponential')

nagen.mod <- list(nag.wei, nagd.wei, nagl.wei, nagdl.wei,
                  nag.exp, nagd.exp, nagl.exp, nagdl.exp)

# just europe
# species
er.wei <- survreg(formula = er.surv ~ 1, dist = 'weibull')
erd.wei <- survreg(formula = er.surv ~ diet, data = er.ecol, dist = 'weibull')
erl.wei <- survreg(formula = er.surv ~ move, data = er.ecol, dist = 'weibull')
erdl.wei <- survreg(formula = er.surv ~ diet + move,
                    data = er.ecol, dist = 'weibull')
er.exp <- survreg(formula = er.surv ~ 1, dist = 'exponential')
erd.exp <- survreg(formula = er.surv ~ diet, data = er.ecol, dist = 'exponential')
erl.exp <- survreg(formula = er.surv ~ move, data = er.ecol, dist = 'exponential')
erdl.exp <- survreg(formula = er.surv ~ diet + move,
                    data = er.ecol, dist = 'exponential')

er.mod <- list(er.wei, erd.wei, erl.wei, erdl.wei,
               er.exp, erd.exp, erl.exp, erdl.exp)
# genus
erg.wei <- survreg(formula = ergen.surv ~ 1, dist = 'weibull')
ergd.wei <- survreg(formula = ergen.surv ~ diet, data = er.genecol, dist = 'weibull')
ergl.wei <- survreg(formula = ergen.surv ~ move, data = er.genecol, dist = 'weibull')
ergdl.wei <- survreg(formula = ergen.surv ~ diet + move,
                     data = er.genecol, dist = 'weibull')
erg.exp <- survreg(formula = ergen.surv ~ 1, dist = 'exponential')
ergd.exp <- survreg(formula = ergen.surv ~ diet, data = er.genecol, dist = 'exponential')
ergl.exp <- survreg(formula = ergen.surv ~ move, data = er.genecol, dist = 'exponential')
ergdl.exp <- survreg(formula = ergen.surv ~ diet + move,
                     data = er.genecol, dist = 'exponential')

ergen.mod <- list(erg.wei, ergd.wei, ergl.wei, ergdl.wei,
                  erg.exp, ergd.exp, ergl.exp, ergdl.exp)
