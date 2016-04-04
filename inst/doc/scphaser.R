## ------------------------------------------------------------------------
##Extract data-structures from the humantcell list
library('scphaser')
invisible(humantcell)
featdata = humantcell[['featdata']]
refcount = humantcell[['refcount']]
altcount = humantcell[['altcount']]
phenodata = humantcell[['phenodata']]

## ------------------------------------------------------------------------
acset = new_acset(featdata, refcount, altcount, phenodata)
lapply(acset, dim)

## ------------------------------------------------------------------------
nminvar = 2
acset = filter_feat_nminvar(acset, nminvar)
lapply(acset, dim)

## ------------------------------------------------------------------------
min_acount = 3
fc = 3
acset = call_gt(acset, min_acount, fc)
lapply(acset, dim)

## ------------------------------------------------------------------------
acset_rnd = racset(acset, type = 'gt')

## ------------------------------------------------------------------------

## ------------------------------------------------------------------------
nmincells = 3
nminvar = 2
acset = filter_acset(acset, nmincells, nminvar)
lapply(acset, dim)

## ------------------------------------------------------------------------
length(unique(acset$featdata$feat))

## ------------------------------------------------------------------------
##acset = phase(acset, input = 'ac', weigh = FALSE, method = 'exhaust')
##acset = phase(acset, input = 'ac', weigh = FALSE, method = 'cluster')
##acset = phase(acset, input = 'ac', weigh = TRUE, method = 'exhaust')
##acset = phase(acset, input = 'ac', weigh = TRUE, method = 'cluster')

##acset = phase(acset, input = 'gt', weigh = TRUE, method = 'exhaust')
##acset = phase(acset, input = 'gt', weigh = TRUE, method = 'cluster')
##acset = phase(acset, input = 'gt', weigh = FALSE, method = 'cluster')
acset = phase(acset, input = 'gt', weigh = FALSE, method = 'exhaust', verbosity = 0)

lapply(acset, dim)

##Number of variants changed
length(acset[['varflip']])

##scores
length(unique(acset[['featdata']][, 'feat']))
length(unique(acset[['score']]))

## ------------------------------------------------------------------------
##gt concordance before and after phasing
acset = set_gt_conc(acset)

##print concordance before and after phasing
acset$gt_conc$conc$feat2ncell
acset$gt_phased_conc$conc$feat2ncell

##print inconcordance
acset$gt_conc$notconc$feat2ncell
acset$gt_phased_conc$notconc$feat2ncell

## ------------------------------------------------------------------------

##filter
nmincells = 3
nminvar = 2
acset_rnd = filter_acset(acset_rnd, nmincells, nminvar)
lapply(acset_rnd, dim)
length(unique(acset_rnd$featdata$feat))

##phase
acset_rnd = phase(acset_rnd, input = 'gt', weigh = FALSE, method = 'exhaust', verbosity = 0)

##gt concordance before and after phasing
acset_rnd = set_gt_conc(acset_rnd)    

##print concordance before and after phasing
acset_rnd$gt_conc$conc$feat2ncell
acset_rnd$gt_phased_conc$conc$feat2ncell

##print inconcordance after phasing
acset_rnd$gt_conc$notconc$feat2ncell
acset_rnd$gt_phased_conc$notconc$feat2ncell

## ----label, fig.show='hold',fig.cap='Genotype concordance between two variants within each gene.'----
plot_conc(acset)
plot_conc(acset_rnd)

