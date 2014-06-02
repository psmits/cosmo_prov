library(survival)
library(MuMIn)

source('../R/na_surv.r')
source('../R/eur_surv.r')

# just north america
# species
na.km <- survfit(formula = na.surv ~ 1)
na.kmd <- survfit(formula = na.surv ~ diet, data = na.ecol)
na.kml <- survfit(formula = na.surv ~ move, data = na.ecol)
na.kmdl <- survfit(formula = na.surv ~ diet + move, data = na.ecol)

# genus
nag.km <- survfit(formula = nagen.surv ~ 1)
nag.kmd <- survfit(formula = nagen.surv ~ diet, data = na.genecol)
nag.kml <- survfit(formula = nagen.surv ~ move, data = na.genecol)
nag.kmdl <- survfit(formula = nagen.surv ~ diet + move, data = na.genecol)

# europe
# species
er.km <- survfit(formula = er.surv ~ 1)
er.kmd <- survfit(formula = er.surv ~ diet, data = er.ecol)
er.kml <- survfit(formula = er.surv ~ move, data = er.ecol)
er.kmdl <- survfit(formula = er.surv ~ diet + move, data = er.ecol)

# genus
erg.km <- survfit(formula = ergen.surv ~ 1)
erg.kmd <- survfit(formula = ergen.surv ~ diet, data = er.genecol)
erg.kml <- survfit(formula = ergen.surv ~ move, data = er.genecol)
erg.kmdl <- survfit(formula = ergen.surv ~ diet + move, data = er.genecol)
