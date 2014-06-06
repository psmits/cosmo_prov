library(survival)
library(MuMIn)
library(xtable)

source('../R/make_table.r')
source('../R/model_sel.r')

source('../R/surv_parametric.r')

# model selection tables
# species
na.table <- surv.tab(na.mod, 'tab:na')
print.xtable(na.table, file = '../doc/namod_tabs.tex',
             hline.after = 0,
             include.rownames = FALSE)

er.table <- surv.tab(er.mod, 'tab:er')
print.xtable(er.table, file = '../doc/ermod_tabs.tex',
             hline.after = 0,
             include.rownames = FALSE)

# genera
nagen.table <- surv.tab(nagen.mod, 'tab:nag')
print.xtable(nagen.table, file = '../doc/nagenmod_tabs.tex',
             hline.after = 0,
             include.rownames = FALSE)

ergen.table <- surv.tab(ergen.mod, 'tab:erg')
print.xtable(ergen.table, file = '../doc/ergenmod_tabs.tex',
             hline.after = 0,
             include.rownames = FALSE)

# relative variable importance
na.imp <- var.imp(na.mod)
er.imp <- var.imp(er.mod)

na.med <- ddply(Reduce(rbind, pna.mod), .(pred), summarize, 
                baseline = median(imp))
er.med <- ddply(Reduce(rbind, per.mod), .(pred), summarize, 
                baseline = median(imp))

nagen.imp <- var.imp(nagen.mod)
ergen.imp <- var.imp(ergen.mod)

nagen.med <- ddply(Reduce(rbind, pnagen.mod), .(pred), summarize, 
                   baseline = median(imp))
ergen.med <- ddply(Reduce(rbind, pergen.mod), .(pred), summarize, 
                   baseline = median(imp))

wei <- laply(na.mod, function(x) x$dist == 'weibull')
na.shape <- lapply(na.mod[wei], ext.shape)
wei <- laply(er.mod, function(x) x$dist == 'weibull')
er.shape <- lapply(er.mod[wei], ext.shape)

wei <- laply(nagen.mod, function(x) x$dist == 'weibull')
nagen.shape <- lapply(nagen.mod[wei], ext.shape)
wei <- laply(ergen.mod, function(x) x$dist == 'weibull')
ergen.shape <- lapply(ergen.mod[wei], ext.shape)
