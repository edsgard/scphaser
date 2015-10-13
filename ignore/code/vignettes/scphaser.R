

## /*
##******************
## */
##' ##Create data-structure
##' Create an "acset" data-structure. An acset is a list that at a minimum contains four data-structures:
##'
##' - featdata: A data-frame which maps variant names to features. It must contain two columns named "feat" and "var", and the rownames must be set to the "var" column.
##' - phenodata: A data-frame which annotates the samples. It must contain a column named "sample". If not provided it is created by the "new_acset" function.
##' - refcount: A matrix with allelic counts for the reference allele. The rownames have to match the values in the "var" column in the featdata, and the colnames the values in the phenodata "sample" column.
##' - altcount: A matrix with allelic counts for the alternative allele, with the same row- and col-names as refcount.
##'
##' As part of this package a dataset with synthetic allele counts are provided, called "fakemouse". The dataset is a list containing the four required data-structures to create an "acset".
##+
##Extract data-structures from the fakemouse list
library('scphaser')
invisible(fakemouse)
featdata = fakemouse[['featdata']]
refcount = fakemouse[['refcount']]
altcount = fakemouse[['altcount']]
phenodata = fakemouse[['phenodata']]

##' Create an "acset" data-structure
##+
acset = new_acset(featdata, refcount, altcount, phenodata)
lapply(acset, dim)

##' Randomize original counts, before any filtering is done. This swaps half of the elements btw the ref and alt matrixes. The randomized dataset will be used below.
##+
acset_rnd = racset(acset)
lapply(acset, dim)


## /*
##******************
## */
##' ## Filter on at least n variants per feature
##+
nminvar = 2
acset = filter_nminvar(acset, nminvar)
lapply(acset, dim)


## /*
##******************
## */
##' ##Call genotypes
##+
min_acount = 3
fc = 50
acset = call_gt(acset, min_acount, fc)
lapply(acset, dim)


## /*
##******************
## */
##' ##Filter variants with no presence of both mono-allelic call
##+
acset = filter_nminmono(acset)
lapply(acset, dim)


## /*
##******************
## */
##' ##Filter on at least n variants per feature
##+
nminvar = 2
acset = filter_nminvar(acset, nminvar)
lapply(acset, dim)


## /*
##******************
## */
##' ##Phase
##+
acset = phase(acset)
lapply(acset, dim)

##Number of variants unchanged (FALSE) and changed (TRUE)
table(acset[['varflip']])


## /*
##******************
## */
##' ##Assess genotype concordance
##+
##gt concordance before and after phasing
acset = set_gt_conc(acset)

##print concordance before and after phasing
acset$gt_conc$conc$feat2ncell
acset$gt_phased_conc$conc$feat2ncell

##print inconcordance
acset$gt_conc$notconc$feat2ncell
acset$gt_phased_conc$notconc$feat2ncell


## /*
##******************
## */
##' ## Phasing of a randomized genotype matrix
##+
acset_rnd = filter_nminvar(acset_rnd, nminvar)
acset_rnd = call_gt(acset_rnd, min_acount, fc)
acset_rnd = filter_nminmono(acset_rnd)
acset_rnd = filter_nminvar(acset_rnd, nminvar)
acset_rnd = phase(acset_rnd)

##gt concordance before and after phasing
acset_rnd = set_gt_conc(acset_rnd)    

##print concordance before and after phasing
acset_rnd$gt_conc$conc$feat2ncell
acset_rnd$gt_phased_conc$conc$feat2ncell

##print inconcordance after phasing
acset_rnd$gt_conc$notconc$feat2ncell
acset_rnd$gt_phased_conc$notconc$feat2ncell


## /*
##******************
## */
##' ##Plot
##+ label, fig.show='hold',fig.cap='Genotype concordance between two variants within each gene.'
plot_conc(acset)
plot_conc(acset_rnd)

