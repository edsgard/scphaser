
library('scphaser')
context('Phasing of genotype matrices')

test_that('phase_cluster phases correctly', {

    ##*###
    ##Half of the cells maternal, half paternal
    ##*###
    
    ##create gt matrix
    ncells = 10
    paternal = c(0, 2, 0, 0, 2)
    maternal = c(2, 0, 2, 2, 0)

    ##expected variants to be flipped: (TBD: possibly c(1, 3, 4))
    vars2flip_expected = as.integer(c(2, 5))

    gt = as.matrix(as.data.frame(rep(list(paternal, maternal), ncells / 2)))
    vars = 1:nrow(gt)
    colnames(gt) = 1:ncells
    rownames(gt) = vars

    ##featdata
    nvars = nrow(gt)
    featdata = as.data.frame(matrix(cbind(rep('jfeat', nvars), 1:nvars), ncol = 2, dimnames = list(vars, c('feat', 'var'))))
    
    ##create acset
    acset = new_acset(featdata, gt = gt)

    ##phase
    vars2flip = phase_cluster(acset, vars)

    ##test
    expect_identical(vars2flip, vars2flip_expected)

    
    ##*###
    ##All cells of one and the same phase
    ##*###
    ncells = 10
    paternal = c(0, 2, 0, 0, 2)
    gt = as.matrix(as.data.frame(rep(list(paternal), ncells)))
    
    ##expected variants to be flipped: (TBD: possibly c(1, 3, 4))
    vars2flip_expected = as.integer(c(2, 5))
    
    vars = 1:nrow(gt)
    colnames(gt) = 1:ncells
    rownames(gt) = vars

    ##featdata
    nvars = nrow(gt)
    featdata = as.data.frame(matrix(cbind(rep('jfeat', nvars), 1:nvars), ncol = 2, dimnames = list(vars, c('feat', 'var'))))
    
    ##create acset
    acset = new_acset(featdata, gt = gt)

    ##phase
    vars2flip = phase_cluster(acset, vars)

    ##test
    expect_identical(vars2flip, vars2flip_expected)
    

    ##*###
    ##Pre-phased (matrix with identical values)
    ##*###

    ncells = 10
    paternal = c(0, 0, 0, 0, 0)
    gt = as.matrix(as.data.frame(rep(list(paternal), ncells)))
    
    ##expected variants to be flipped: (TBD: possibly c(1, 3, 4))
    vars2flip_expected = character(0)
    
    vars = 1:nrow(gt)
    colnames(gt) = 1:ncells
    rownames(gt) = vars

    ##featdata
    nvars = nrow(gt)
    featdata = as.data.frame(matrix(cbind(rep('jfeat', nvars), 1:nvars), ncol = 2, dimnames = list(vars, c('feat', 'var'))))
    
    ##create acset
    acset = new_acset(featdata, gt = gt)

    ##phase
    vars2flip = phase_cluster(acset, vars)

    ##test
    expect_identical(vars2flip, vars2flip_expected)
    
    
    ##*###
    ##TBD: test if works with added "noise" in the form of 1's and NAs
    ##*###
    
})

test_that('phase_exhaustive phases correctly', {

    ##*###
    ##Half of the cells maternal, half paternal
    ##*###
    
    ##create gt matrix
    ncells = 10
    paternal = c(0, 2, 0, 0, 2)
    maternal = c(2, 0, 2, 2, 0)

    ##expected variants to be flipped
    vars2flip_expected = as.integer(c(1, 3, 4))

    gt = as.matrix(as.data.frame(rep(list(paternal, maternal), ncells / 2)))
    vars = 1:nrow(gt)
    colnames(gt) = 1:ncells
    rownames(gt) = vars

    ##featdata
    nvars = nrow(gt)
    featdata = as.data.frame(matrix(cbind(rep('jfeat', nvars), 1:nvars), ncol = 2, dimnames = list(vars, c('feat', 'var'))))
    
    ##create acset
    acset = new_acset(featdata, gt = gt)

    ##phase
    vars2flip = phase_exhaustive(acset, vars)

    ##test
    expect_identical(vars2flip, vars2flip_expected)

    
    ##*###
    ##All cells of one and the same phase
    ##*###
    ncells = 10
    paternal = c(0, 2, 0, 0, 2)
    gt = as.matrix(as.data.frame(rep(list(paternal), ncells)))
    
    ##expected variants to be flipped
    vars2flip_expected = as.integer(c(1, 3, 4))
    
    vars = 1:nrow(gt)
    colnames(gt) = 1:ncells
    rownames(gt) = vars

    ##featdata
    nvars = nrow(gt)
    featdata = as.data.frame(matrix(cbind(rep('jfeat', nvars), 1:nvars), ncol = 2, dimnames = list(vars, c('feat', 'var'))))
    
    ##create acset
    acset = new_acset(featdata, gt = gt)

    ##phase
    vars2flip = phase_exhaustive(acset, vars)

    ##test
    expect_identical(vars2flip, vars2flip_expected)
    

    ##*###
    ##Pre-phased (matrix with identical values)
    ##*###

    ncells = 10
    paternal = c(0, 0, 0, 0, 0)
    gt = as.matrix(as.data.frame(rep(list(paternal), ncells)))
    
    ##expected variants to be flipped: (TBD: possibly c(1, 3, 4))
    vars2flip_expected = integer(0)
    
    vars = 1:nrow(gt)
    colnames(gt) = 1:ncells
    rownames(gt) = vars

    ##featdata
    nvars = nrow(gt)
    featdata = as.data.frame(matrix(cbind(rep('jfeat', nvars), 1:nvars), ncol = 2, dimnames = list(vars, c('feat', 'var'))))
    
    ##create acset
    acset = new_acset(featdata, gt = gt)

    ##phase
    vars2flip = phase_exhaustive(acset, vars)

    ##test
    expect_identical(vars2flip, vars2flip_expected)
    
    
    ##*###
    ##TBD: test if works with added "noise" in the form of 1's and NAs
    ##*###
    
})
