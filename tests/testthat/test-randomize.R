
library('scphaser')
context('Randomization of acset data-structure')

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

test_that('racset permutes allele counts with fixed rowmargins', {

    ##create alt and ref matrixes with no shared values
    nels = 100
    nrows = 5
    ncols = nels / nrows
    refcount = matrix(sample(1:nels, nels), nrow = nrows, dimnames = list(1:nrows, 1:ncols))
    altcount = matrix(sample((nels+1):(nels*2), nels), nrow = nrows, dimnames = list(1:nrows, 1:ncols))

    acset = new_acset(rownames(refcount), refcount, altcount)

    ##randomize
    rac = racset(acset, type = 'ac', fixedrowmargin = TRUE)
    
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
    ##check that the rowsums are fixed
    ##*##

    ##expected
    preref_rowsums = apply(refcount, 1, sum)
    prealt_rowsums = apply(altcount, 1, sum)
    
    ##observed
    postref_rowsums = apply(rac[['refcount']], 1, sum)
    postalt_rowsums = apply(rac[['altcount']], 1, sum)
    
    ##test
    expect_identical(postref_rowsums, preref_rowsums)
    expect_identical(postalt_rowsums, prealt_rowsums)
    
})

test_that('racset permutes gt with fixed rowmargins', {

    ##create alt and ref matrixes with no shared values
    acset = gt_matpat()$acset
    gt = acset$gt
    
    ##randomize
    rac = racset(acset, type = 'gt', fixedrowmargin = TRUE)
    gt_rnd = rac$gt
    
    ##*###
    ##check that total number of elements hasn't changed
    ##*###
    
    ##expected
    origsize = length(gt)

    ##observed
    expsize = length(gt_rnd)

    ##test
    expect_identical(expsize, origsize)

    
    ##*##
    ##check that the rowsums are fixed
    ##*##

    ##expected
    pregt_rowsums = apply(gt, 1, table)
    
    ##observed
    postgt_rowsums = apply(gt_rnd, 1, table)
    
    ##test
    expect_identical(postgt_rowsums, pregt_rowsums)
    
})

test_that('racset permutes gt with fixed table sum', {

    ##create alt and ref matrixes with no shared values
    acset = gt_matpat()$acset
    gt = acset$gt
    
    ##randomize
    rac = racset(acset, type = 'gt', fixedrowmargin = FALSE)
    gt_rnd = rac$gt
    
    ##*###
    ##check that total number of elements hasn't changed
    ##*###
    
    ##expected
    origsize = length(gt)

    ##observed
    expsize = length(gt_rnd)

    ##test
    expect_identical(expsize, origsize)

    
    ##*##
    ##check that the rowsums are fixed
    ##*##

    ##expected
    pregt_sum = sum(gt)
    
    ##observed
    postgt_sum = sum(gt_rnd)
    
    ##test
    expect_identical(postgt_sum, pregt_sum)
    
})

test_that('racset permutes allele counts with fixed table margin', {

    ##create alt and ref matrixes with no shared values
    acset = ac_matpat()$acset

    ##randomize
    rac = racset(acset, type = 'ac', fixedrowmargin = FALSE)
    
    ##*###
    ##check that total number of elements hasn't changed
    ##*###
    
    ##expected
    refcount = acset$refcount
    altcount = acset$altcount
    origsize = length(refcount)

    ##observed
    refsize = length(rac[['refcount']])
    altsize = length(rac[['altcount']])

    ##test
    expect_identical(refsize, origsize)
    expect_identical(altsize, origsize)

    
    ##*##
    ##check that the rowsums are not fixed
    ##*##

    ##expected
    preref_rowsums = apply(refcount, 1, sum)
    prealt_rowsums = apply(altcount, 1, sum)
    
    ##observed
    postref_rowsums = apply(rac[['refcount']], 1, sum)
    postalt_rowsums = apply(rac[['altcount']], 1, sum)
    
    ##test
    expect_false(identical(postref_rowsums, preref_rowsums))
    expect_false(identical(postalt_rowsums, prealt_rowsums))

    ##*##
    ##check that table sum is fixed
    ##*##

    ##expected    
    preref_sum = sum(refcount)
    prealt_sum = sum(altcount)
    
    ##observed
    postref_sum = sum(rac[['refcount']])
    postalt_sum = sum(rac[['altcount']])
    
    ##test
    expect_identical(postref_sum, preref_sum)
    expect_identical(postalt_sum, prealt_sum)
})


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
    verbosity = 0
    
    nperm = 100

    ##preproc
    acset = filter_nminvar(acset, nminvar)
    
    ##phase
    acset = phase(acset, input, weigh, method, nvars_max, verbosity)
    
    ##get p-value
    pval = get_phase_pval(acset, nperm)

    
})
