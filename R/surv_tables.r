library(survival)
library(MuMIn)
library(xtable)

source('../R/make_table.r')

source('../R/surv_analysis.r')

na.table <- surv.tab(na.mod, 'tab:na')

print.xtable(na.table, file = '../doc/namod_tabs.tex',
             hline.after = 0,
             include.rownames = FALSE)

er.table <- surv.tab(er.mod, 'tab:er')

print.xtable(na.table, file = '../doc/ermod_tabs.tex',
             hline.after = 0,
             include.rownames = FALSE)
