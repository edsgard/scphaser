
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

gt_prephased <- function(){
    
    ##*###
    ##Pre-phased (matrix with identical values)
    ##*###

    ncells = 10
    paternal = c(0, 2, 0, 2)
    gt = as.matrix(as.data.frame(rep(list(paternal), ncells)))
        
    vars = 1:nrow(gt)
    colnames(gt) = 1:ncells
    rownames(gt) = vars

    ##featdata
    nvars = nrow(gt)
    featdata = as.data.frame(matrix(cbind(rep('jfeat', nvars), as.character(1:nvars)), ncol = 2, dimnames = list(vars, c('feat', 'var'))))
    
    ##create acset
    acset = new_acset(featdata, gt = gt)

    return(acset = acset)
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
    verbosity = 0
    
    nperm = 100

    ##preproc
    acset = filter_feat_nminvar(acset, nminvar)
    
    ##phase
    acset = phase(acset, input, weigh, method, nvars_max, verbosity)
    
    ##get p-value
    pval = get_phase_pval(acset, nperm)


    ##*###
    ##Pre-phased
    ##*###
    acset = gt_prephased()

    ##params
    nminvar = 2
    input = 'gt'
    weigh = FALSE
    method = 'exhaust'
    nvars_max = 10
    verbosity = 0
    
    nperm = 100

    ##preproc
    acset = filter_feat_nminvar(acset, nminvar)
    
    ##phase
    acset = phase(acset, input, weigh, method, nvars_max, verbosity)
    
    ##get p-value
    pval = get_phase_pval(acset, nperm, fixedrowmargin = TRUE)
    
    pval = get_phase_pval(acset, nperm, fixedrowmargin = FALSE)
    
})

gt_matpat_varscore <- function(){

    ##create gt matrix
    ncells = 8
    paternal = c(0, 2, 1, 0, 0, 2)
    maternal = c(2, 0, 1, 2, 2, 0)
    
    gt = as.matrix(as.data.frame(rep(list(paternal, maternal), ncells / 2)))

    ##add row with only NAs
    gt = rbind(gt, rep(NA, ncol(gt)))

    ##add row with only 1s
    gt = rbind(gt, rep(1, ncol(gt)))

    ##add column with only NAs
    gt = cbind(gt, rep(NA, nrow(gt)))

    ##add column with majority NAs, where 2 should become the reference
    gt = cbind(gt, c(rep(NA, 3), rep(1, 2), 0, rep(2, 2)))
        
    vars = 1:nrow(gt)
    colnames(gt) = 1:ncol(gt)
    rownames(gt) = vars

    ##featdata
    nvars = nrow(gt)
    featdata = as.data.frame(matrix(cbind(rep('jfeat', nvars), 1:nvars), ncol = 2, dimnames = list(vars, c('feat', 'var'))))
    
    ##create acset
    acset = new_acset(featdata, gt = gt)

    ##expected
    expect = matrix(c(0.8, 0.2, 0,
        0.8, 0.2, 0,
        0, 0.2, 0.8,
        0.8, 0.1, 0.1,
        0.8, 0.1, 0.1,
        0.9, 0.1, 0,
        0.1, 0.9, 0,
        0.1, 0.1, 0.8), ncol = 3, byrow = T)
    colnames(expect) = c('gt', 'na', 'bi')
    rownames(expect) = 1:nrow(expect)

    return(list(acset = acset, expect = expect))
}

test_that('variant score calculation works', {

    ##dataset
    acset_exp = gt_matpat_varscore()
    acset = acset_exp$acset
    expect = acset_exp$expect
    
    ##params
    nminvar = 2
    input = 'gt'
    weigh = FALSE
    method = 'exhaust'
    nvars_max = 10
    verbosity = 0
    
    ##preproc
    acset = filter_feat_nminvar(acset, nminvar)
    
    ##phase
    acset = phase(acset, input, weigh, method, nvars_max, verbosity)

    ##variant scores
    varscore = score_var(acset)

    ##test
    expect_identical(varscore, expect)
})

