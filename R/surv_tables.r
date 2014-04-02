library(survival)
library(MuMIn)
library(xtable)

source('../R/make_table.r')

source('../R/surv_parametric.r')

na.table <- surv.tab(na.mod, 'tab:na')

print.xtable(na.table, file = '../doc/namod_tabs.tex',
             hline.after = 0,
             include.rownames = FALSE)

best.table <- surv.tab(best.mod, 'tab:best')

er.table <- surv.tab(er.mod, 'tab:er')

print.xtable(er.table, file = '../doc/ermod_tabs.tex',
             hline.after = 0,
             include.rownames = FALSE)


nagen.table <- surv.tab(nagen.mod, 'tab:nag')

print.xtable(nagen.table, file = '../doc/nagenmod_tabs.tex',
             hline.after = 0,
             include.rownames = FALSE)

ergen.table <- surv.tab(ergen.mod, 'tab:erg')

print.xtable(ergen.table, file = '../doc/ergenmod_tabs.tex',
             hline.after = 0,
             include.rownames = FALSE)
