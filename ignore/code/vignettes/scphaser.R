

## /*
##******************
## */
##' ##Create data-structure
##' Create an "acset" data-structure. An acset is a list that at a minimum contains four data-structures:
##'
##' - featdata: A data-frame which maps variant names to feature names and which contains reference and alternative alleles. It must contain four columns named "feat", "var", "ref" and "alt", and the rownames must be set to the "var" column.
##' - phenodata: A data-frame which annotates the samples. It must contain a column named "sample". If not provided it is created by the "new_acset" function.
##' - refcount: A matrix with allelic counts for the reference allele. The rownames have to match the values in the "var" column in the featdata, and the colnames the values in the phenodata "sample" column.
##' - altcount: A matrix with allelic counts for the alternative allele, with the same row- and col-names as refcount.
##'
##' As part of this package a dataset with allele counts from human are provided, called "marinov". The dataset is a list containing the four required data-structures to create an "acset".
##+
##Extract data-structures from the marinov list
library('scphaser')
invisible(marinov)
featdata = marinov[['featdata']]
refcount = marinov[['refcount']]
altcount = marinov[['altcount']]
phenodata = marinov[['phenodata']]

##' Create an "acset" data-structure
##+
acset = new_acset(featdata, refcount, altcount, phenodata)

##' Print the elements of the acset data-structure and the dimensions of these, illustrating the number of variants and cells.
##+
lapply(acset, dim)

##' Print the number of genes
##+
length(unique(acset$featdata$feat))


## /*
##******************
## */
##' ## Filter on at least n variants per feature
##' Features with less than "nminvar" variants are removed
##+
nminvar = 2
acset = filter_feat_nminvar(acset, nminvar)

##' Print the element dimensions after filtering out genes with less than 2 variants
lapply(acset, dim)


## /*
##******************
## */
##' ##Call genotypes
##' Transcribed genotypes are called as 2 or 0 if there at least "min_acount" reads and the fold-change >= 3 or <= 1/3, where fold-change = alternative allele count / reference allele count. For entries that do not meet the criteria, such as when bi-allelic expression close to a 50/50 expression ratio between the alleles, the genotype is set to 1.
##+
min_acount = 3
fc = 3
acset = call_gt(acset, min_acount, fc)

##' Print the elements after calling genotypes. Note that a new element called "gt" was created, containing the called genotypes.
lapply(acset, dim)


##' Randomize original counts. The randomized dataset will be used below.
##+
acset_rnd = racset(acset, type = 'gt')


## /*
##******************
## */
##' ## Filter variants
##+

##' Filter variants on having at least "nmincells" cells with monoallelic calls in at least "nminvar" variants within a feature. After variant filtering, filter features on having at least "nminvar" variants.
##+
nmincells = 5
nminvar = 2
acset = filter_acset(acset, nmincells, nminvar)

##' Number of variants after filtering
##+
nrow(acset$featdata)

##' Number of features after filtering
##+
length(unique(acset$featdata$feat))


## /*
##******************
## */
##' ##Phase
##' There are three phasing arguments that can be provided to the phasing function, each with two possible values. Below we call it using genotypes, exhaustive clustering and without weighing and also illustrate calling using the other possible values (first commented lines).
##+
##acset = phase(acset, input = 'ac', weigh = FALSE, method = 'exhaust')
##acset = phase(acset, input = 'ac', weigh = FALSE, method = 'pam')
##acset = phase(acset, input = 'ac', weigh = TRUE, method = 'exhaust')
##acset = phase(acset, input = 'ac', weigh = TRUE, method = 'pam')

##acset = phase(acset, input = 'gt', weigh = TRUE, method = 'exhaust')
##acset = phase(acset, input = 'gt', weigh = TRUE, method = 'pam')
##acset = phase(acset, input = 'gt', weigh = FALSE, method = 'pam')
acset = phase(acset, input = 'gt', weigh = FALSE, method = 'exhaust', verbosity = 0)

##' As an overview of the elements in the acset datastructure we print the dimension of each element
lapply(acset, dim)

##' The haplotype output is contained in an element that was added by the phasing function and is named "phasedfeat"
head(acset[['phasedfeat']])

##' The variants that have been swapped compared to the input reference and alternative alleles are contained in "varflip". The number of swapped variants is the length of this vector.
length(acset[['varflip']])


## /*
##******************
## */
##' ##Assess genotype concordance
##' To assess the success of the phasing one can calculate the degree of variability remaining if all cells with haplotype 2 are set to haplotype 1. As a rough measure of this we here calculate the variability as the number of cells that differ from the inferred haplotype for every gene with two variants. The differing number of cells per gene we here denote inconcordance.
##+
##get gt concordance before and after phasing
acset = set_gt_conc(acset)

##' Inconcordance before phasing
table(acset$gt_conc$notconc$feat2ncell)

##' Inconcordance after phasing
table(acset$gt_phased_conc$notconc$feat2ncell)


## /*
##******************
## */
##' ## Phasing of a randomized genotype matrix

##' Filtering
##+
nmincells = 3
nminvar = 2
acset_rnd = filter_acset(acset_rnd, nmincells, nminvar)

##' Dimensions after filtering
lapply(acset_rnd, dim)

##' Number of genes after filtering
length(unique(acset_rnd$featdata$feat))

##' Phasing
##+
acset_rnd = phase(acset_rnd, input = 'gt', weigh = FALSE, method = 'exhaust', verbosity = 0)

##' Get genotype matrix concordance before and after phasing
acset_rnd = set_gt_conc(acset_rnd)    

##' Inconcordance before phasing
table(acset_rnd$gt_conc$notconc$feat2ncell)

##' Inconcordance after phasing
table(acset_rnd$gt_phased_conc$notconc$feat2ncell)


## /*
##******************
## */
##' ##Plot
##' The removal of variability from the cell distribution after phasing indicates successful phasing. The right subfigure show that phasing cannot be done if the genotype matrix is scrambled, since the transcribed genotype of each cell then do not anylonger correspond to one of two underlying haplotype states.
##+ label, fig.show='hold',fig.cap='Genotype concordance between two variants within each gene.'
plot_conc(acset)
plot_conc(acset_rnd)
