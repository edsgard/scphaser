## ------------------------------------------------------------------------
##Extract data-structures from the marinov list
library('scphaser')
invisible(marinov)
featdata = marinov[['featdata']]
refcount = marinov[['refcount']]
altcount = marinov[['altcount']]
phenodata = marinov[['phenodata']]

## ------------------------------------------------------------------------
acset = new_acset(featdata, refcount, altcount, phenodata)

## ------------------------------------------------------------------------
lapply(acset, dim)

## ------------------------------------------------------------------------
length(unique(acset$featdata$feat))

## ------------------------------------------------------------------------
nminvar = 2
acset = filter_feat_nminvar(acset, nminvar)

## ------------------------------------------------------------------------
lapply(acset, dim)

## ------------------------------------------------------------------------
min_acount = 3
fc = 3
acset = call_gt(acset, min_acount, fc)

## ------------------------------------------------------------------------
lapply(acset, dim)

## ------------------------------------------------------------------------
acset_rnd = racset(acset, type = 'gt')

## ------------------------------------------------------------------------

## ------------------------------------------------------------------------
nmincells = 5
nminvar = 2
acset = filter_acset(acset, nmincells, nminvar)

## ------------------------------------------------------------------------
nrow(acset$featdata)

## ------------------------------------------------------------------------
length(unique(acset$featdata$feat))

## ------------------------------------------------------------------------
##acset = phase(acset, input = 'ac', weigh = FALSE, method = 'exhaust')
##acset = phase(acset, input = 'ac', weigh = FALSE, method = 'pam')
##acset = phase(acset, input = 'ac', weigh = TRUE, method = 'exhaust')
##acset = phase(acset, input = 'ac', weigh = TRUE, method = 'pam')

##acset = phase(acset, input = 'gt', weigh = TRUE, method = 'exhaust')
##acset = phase(acset, input = 'gt', weigh = TRUE, method = 'pam')
##acset = phase(acset, input = 'gt', weigh = FALSE, method = 'pam')
acset = phase(acset, input = 'gt', weigh = FALSE, method = 'exhaust', verbosity = 0)

## ------------------------------------------------------------------------
lapply(acset, dim)

## ------------------------------------------------------------------------
head(acset[['phasedfeat']])

## ------------------------------------------------------------------------
length(acset[['varflip']])

## ------------------------------------------------------------------------
##get gt concordance before and after phasing
acset = set_gt_conc(acset)

## ------------------------------------------------------------------------
table(acset$gt_conc$notconc$feat2ncell)

## ------------------------------------------------------------------------
table(acset$gt_phased_conc$notconc$feat2ncell)

## ------------------------------------------------------------------------
nmincells = 3
nminvar = 2
acset_rnd = filter_acset(acset_rnd, nmincells, nminvar)

## ------------------------------------------------------------------------
lapply(acset_rnd, dim)

## ------------------------------------------------------------------------
length(unique(acset_rnd$featdata$feat))

## ------------------------------------------------------------------------
acset_rnd = phase(acset_rnd, input = 'gt', weigh = FALSE, method = 'exhaust', verbosity = 0)

## ------------------------------------------------------------------------
acset_rnd = set_gt_conc(acset_rnd)    

## ------------------------------------------------------------------------
table(acset_rnd$gt_conc$notconc$feat2ncell)

## ------------------------------------------------------------------------
table(acset_rnd$gt_phased_conc$notconc$feat2ncell)

## ----label, fig.show='hold',fig.cap='Genotype concordance between two variants within each gene.'----
plot_conc(acset)
plot_conc(acset_rnd)

