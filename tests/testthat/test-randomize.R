
library('scphaser')
context('Randomization of acset data-structure')

test_that('racset swaps allele counts', {

    ##create alt and ref matrixes with no shared values
    nels = 100
    nrows = 5
    ncols = nels / nrows
    refcount = matrix(sample(1:nels, nels), nrow = nrows, dimnames = list(1:nrows, 1:ncols))
    altcount = matrix(sample((nels+1):(nels*2), nels), nrow = nrows, dimnames = list(1:nrows, 1:ncols))

    acset = new_acset(1:rownames(refcount), refcount, altcount)

    ##randomize
    rac = racset(acset)

    
    ##*###
    ##check that total number of elements hasn't changed
    ##*###
    
    ##expected
    origsize = length(refcount)

    ##observed
    refsize = length(rac[['refcount']])
    altsize = length(rac[['altcount']])

    ##test
    expect_identical(refsize, origsize)
    expect_identical(altsize, origsize)

    
    ##*##
    ##check that half of the elements are identical to the original
    ##*##

    ##expected
    nswapped = as.integer(floor(origsize / 2))
    nkept = origsize - nswapped

    ##observed
    nref_kept = length(which(rac[['refcount']] == refcount))
    nref_swapped = length(which(rac[['refcount']] == altcount))
    nalt_kept = length(which(rac[['altcount']] == altcount))
    nalt_swapped = length(which(rac[['altcount']] == refcount))

    ##test
    expect_identical(nref_swapped, nswapped)
    expect_identical(nalt_swapped, nswapped)
    expect_identical(nref_kept, nkept)
    expect_identical(nalt_kept, nkept)    
    
})

gt_matpat <- function(){

    ##*###
    ##Half of the cells maternal, half paternal
    ##*###

    ##create gt matrix
    ncells = 10
    paternal = c(0, 2, 0, 0, 2)
    maternal = c(2, 0, 2, 2, 0)

    ##expected variants to be flipped: (TBD: possibly c(1, 3, 4))
    vars2flip_expected = list(as.character(c(2, 5)), as.character(c(1, 3, 4)))

    gt = as.matrix(as.data.frame(rep(list(paternal, maternal), ncells / 2)))
    vars = 1:nrow(gt)
    colnames(gt) = 1:ncells
    rownames(gt) = vars

    ##featdata
    nvars = nrow(gt)
    featdata = as.data.frame(matrix(cbind(rep('jfeat', nvars), 1:nvars), ncol = 2, dimnames = list(vars, c('feat', 'var'))))
    
    ##create acset
    acset = new_acset(featdata, gt = gt)

    return(list(acset = acset, exp = vars2flip_expected))
}

ac_matpat <- function(){

    ##*###
    ##Half of the cells maternal, half paternal
    ##*###

    ##allele count matrixes
    ncells = 10
    ref = c(0, 10, 1, 5, 2, 9)
    alt = c(0, 2, 5, 5, 6, 3)    
    vars2flip_exp = list(as.character(c(3, 5)), as.character(c(2, 6)))
    
    refcount = as.matrix(as.data.frame(rep(list(ref, alt), ncells / 2)))
    altcount = as.matrix(as.data.frame(rep(list(alt, ref), ncells / 2)))
    
    vars = 1:nrow(refcount)
    samples = 1:ncells
    colnames(refcount) = samples
    rownames(refcount) = vars
    colnames(altcount) = samples
    rownames(altcount) = vars
    
    ##featdata
    nvars = length(vars)
    featdata = as.data.frame(matrix(cbind(rep('jfeat', nvars), 1:nvars), ncol = 2, dimnames = list(vars, c('feat', 'var'))), stringsAsFactors = FALSE)    
    
    ##create acset
    acset = new_acset(featdata, altcount = altcount, refcount = refcount)

    ##call genotype
    min_acount = 3
    fc = 3 #To avoid for example refmap biases. fc == 3 corresponds to 75/25 ratio. Some of the first with expression of both alleles are then: 3/1, 6/2, 9/3.
    acset = call_gt(acset, min_acount, fc)
    lapply(acset, dim)

    ##set weights
    acset = set_aseweights(acset)
    
    return(list(acset = acset, exp = vars2flip_exp))
}

test_that('p-value calculation works', {

    ##dataset
    acset_exp = ac_matpat()
    acset = acset_exp$acset

    ##params
    nminvar = 2
    input = 'gt'
    weigh = FALSE
    method = 'exhaust'
    nvars_max = 10
    
    nperm = 100

    ##preproc
    acset = filter_nminvar(acset, nminvar)
    
    ##phase
    acset = phase(acset, input, weigh, method, nvars_max)
    
    ##get p-value
    pval = get_phase_pval(acset, nperm)

    
})
