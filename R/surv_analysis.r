library(survival)
library(MuMIn)
library(xtable)

source('../R/na_surv.r')
source('../R/eur_surv.r')

aic.wts <- function(aic) {
  dels <- aic - min(aic)
  rel <- exp(-0.5 * dels)

  rel / sum(rel)
}

# just north america
na.km <- survfit(formula = na.surv ~ 1)
na.kmd <- survfit(formula = na.surv ~ diet, data = na.ecol)
na.kml <- survfit(formula = na.surv ~ move, data = na.ecol)
na.kmdl <- survfit(formula = na.surv ~ diet + move, data = na.ecol)

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

# just europe
er.km <- survfit(formula = er.surv ~ 1)
er.kmd <- survfit(formula = er.surv ~ diet, data = er.ecol)
er.kml <- survfit(formula = er.surv ~ move, data = er.ecol)
er.kmdl <- survfit(formula = er.surv ~ diet + move, data = er.ecol)

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

# combined
