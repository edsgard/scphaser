
library('scphaser')
context('Phasing of genotype matrices')

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
    featdata = as.data.frame(matrix(cbind(rep('jfeat', nvars), 1:nvars), ncol = 2, dimnames = list(vars, c('feat', 'var'))))
    
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

gt_samephase <- function(){
    
    ##*###
    ##All cells of one and the same phase
    ##*###
    ncells = 10
    paternal = c(0, 2, 0, 0, 2)
    gt = as.matrix(as.data.frame(rep(list(paternal), ncells)))
    
    ##expected variants to be flipped: (TBD: possibly c(1, 3, 4))
    vars2flip_exp = list(as.character(c(2, 5)), as.character(c(1, 3, 4)))
    
    vars = 1:nrow(gt)
    colnames(gt) = 1:ncells
    rownames(gt) = vars

    ##featdata
    nvars = nrow(gt)
    featdata = as.data.frame(matrix(cbind(rep('jfeat', nvars), 1:nvars), ncol = 2, dimnames = list(vars, c('feat', 'var'))))
    
    ##create acset
    acset = new_acset(featdata, gt = gt)

    return(list(acset = acset, exp = vars2flip_exp))    
}

ac_samephase <- function(){

    ##*###
    ##All cells of one and the same phase
    ##*###
    
    ##allele count matrixes
    ncells = 10
    ref = c(0, 10, 1, 5, 2, 3)
    alt = c(0, 3, 1, 5, 2, 9)    
    vars2flip_exp = list(as.character(2), as.character(6))
    
    refcount = as.matrix(as.data.frame(rep(list(ref), ncells)))
    altcount = as.matrix(as.data.frame(rep(list(alt), ncells)))
    
    vars = 1:nrow(refcount)
    samples = 1:ncells
    colnames(refcount) = samples
    rownames(refcount) = vars
    colnames(altcount) = samples
    rownames(altcount) = vars
    
    ##featdata
    nvars = length(vars)
    featdata = as.data.frame(matrix(cbind(rep('jfeat', nvars), 1:nvars), ncol = 2, dimnames = list(vars, c('feat', 'var'))))
    
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

gt_prephased <- function(){
    
    ##*###
    ##Pre-phased (matrix with identical values)
    ##*###

    ncells = 10
    paternal = c(0, 0, 0, 0, 0)
    gt = as.matrix(as.data.frame(rep(list(paternal), ncells)))
    
    ##expected variants to be flipped: (TBD: possibly c(1, 3, 4))
    vars2flip_exp = character(0)
    
    vars = 1:nrow(gt)
    colnames(gt) = 1:ncells
    rownames(gt) = vars

    ##featdata
    nvars = nrow(gt)
    featdata = as.data.frame(matrix(cbind(rep('jfeat', nvars), as.character(1:nvars)), ncol = 2, dimnames = list(vars, c('feat', 'var'))))
    
    ##create acset
    acset = new_acset(featdata, gt = gt)

    return(list(acset = acset, exp = vars2flip_exp))
}

ac_prephased <- function(){
    
    ##*###
    ##Pre-phased (matrix with identical values)
    ##*###

    ##allele count matrixes
    ncells = 10
    ref = c(0, 10, 1, 5, 2, 9)
    alt = c(0, 3, 1, 5, 2, 3)    
    vars2flip_exp = list(character(0), as.character(c(2, 6)))
    
    refcount = as.matrix(as.data.frame(rep(list(ref), ncells)))
    altcount = as.matrix(as.data.frame(rep(list(alt), ncells)))
    
    vars = 1:nrow(refcount)
    samples = 1:ncells
    colnames(refcount) = samples
    rownames(refcount) = vars
    colnames(altcount) = samples
    rownames(altcount) = vars
    
    ##featdata
    nvars = length(vars)
    featdata = as.data.frame(matrix(cbind(rep('jfeat', nvars), as.character(1:nvars)), ncol = 2, dimnames = list(vars, c('feat', 'var'))))
    
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

ac_weighedtoflip <- function(){
    
    ##Weigh the difference heavy, such that a flip should take place
    
    ##allele count matrixes
    refcount = matrix(c(0, 100, 5, 10, 3, 5), nrow = 3)
    altcount = matrix(c(10, 0, 5, 0, 0, 5), nrow = 3)
    vars2flip_exp = list(as.character(1), as.character(2))
    
    ncells = ncol(refcount)
    vars = 1:nrow(refcount)
    samples = 1:ncells
    colnames(refcount) = samples
    rownames(refcount) = vars
    colnames(altcount) = samples
    rownames(altcount) = vars
    
    ##featdata
    nvars = length(vars)
    featdata = as.data.frame(matrix(cbind(rep('jfeat', nvars), as.character(1:nvars)), ncol = 2, dimnames = list(vars, c('feat', 'var'))))
    
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

ac_weighedtonotflip <- function(){

    ##Weigh the equal gt heavy, such that no flip should take place (or all is flipped)
    
    ##allele count matrixes
    refcount = matrix(c(0, 3, 5, 10, 100, 5), nrow = 3)
    altcount = matrix(c(10, 0, 5, 0, 0, 5), nrow = 3)
    vars2flip_exp = list(as.character(), as.character(c(1, 2)))
    
    ncells = ncol(refcount)
    vars = 1:nrow(refcount)
    samples = 1:ncells
    colnames(refcount) = samples
    rownames(refcount) = vars
    colnames(altcount) = samples
    rownames(altcount) = vars
    
    ##featdata
    nvars = length(vars)
    featdata = as.data.frame(matrix(cbind(rep('jfeat', nvars), as.character(1:nvars)), ncol = 2, dimnames = list(vars, c('feat', 'var'))))
    
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

test_that('can phase gt where half cells maternal and half paternal', {

    ##create gt matrix, providing gt as input
    acset_exp = gt_matpat()
    acset = acset_exp$acset
    exp_out = acset_exp$exp
    vars = rownames(acset[['gt']])
    
    ##exhaust_gt
    vars2flip = phase_exhaust_gt(acset, vars)
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))

    ##cluster_gt
    vars2flip = phase_cluster_gt(acset, vars)
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))

    ##create acset with ac counts as input
    acset_exp = ac_matpat()
    acset = acset_exp$acset
    exp_out = acset_exp$exp
    vars = rownames(acset[['gt']])
        
    ##wexhaust_gt
    vars2flip = wphase_exhaust_gt(acset, vars)
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))

    ##wcluster_gt
    vars2flip = wphase_cluster_gt(acset, vars)
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))            
})

test_that('can phase ase where half cells maternal and half paternal', {

    ##create acset with ac counts as input
    acset_exp = ac_matpat()
    acset = acset_exp$acset
    exp_out = acset_exp$exp
    vars = rownames(acset[['gt']])

    ##exhaust_ac
    vars2flip = phase_exhaust_ase(acset, vars)
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))

    ##wexhaust_ac
    vars2flip = wphase_exhaust_ase(acset, vars)
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))
    
    ##cluster_ac
    vars2flip = phase_cluster_ase(acset, vars)
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))
    
    ##wcluster_ac
    vars2flip = wphase_cluster_ase(acset, vars)
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))    
})

test_that('can phase gt where cells have the same phase', {

    ##create gt matrix, providing gt as input
    acset_exp = gt_samephase()
    acset = acset_exp$acset
    exp_out = acset_exp$exp
    vars = rownames(acset[['gt']])
    
    ##exhaust_gt
    vars2flip = phase_exhaust_gt(acset, vars)
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))

    ##cluster_gt
    vars2flip = phase_cluster_gt(acset, vars)
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))

    ##create acset with ac counts as input
    acset_exp = ac_samephase()
    acset = acset_exp$acset
    exp_out = acset_exp$exp
    vars = rownames(acset[['gt']])
        
    ##wexhaust_gt
    vars2flip = wphase_exhaust_gt(acset, vars)
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))

    ##wcluster_gt
    vars2flip = wphase_cluster_gt(acset, vars)
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))
})

test_that('can phase ase where cells have the same phase', {

    ##create acset with ac counts as input
    acset_exp = ac_samephase()
    acset = acset_exp$acset
    exp_out = acset_exp$exp
    vars = rownames(acset[['gt']])

    ##exhaust_ac
    vars2flip = phase_exhaust_ase(acset, vars)
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))

    ##wexhaust_ac
    vars2flip = wphase_exhaust_ase(acset, vars)
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))
    
    ##cluster_ac
    vars2flip = phase_cluster_ase(acset, vars)
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))
    
    ##wcluster_ac
    vars2flip = wphase_cluster_ase(acset, vars)
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))    
})

test_that('all or nothing change when gt have been prephased', {

    ##create gt matrix, providing gt as input
    acset_exp = gt_prephased()
    acset = acset_exp$acset
    exp_out = acset_exp$exp
    vars = rownames(acset[['gt']])
    
    ##exhaust_gt
    vars2flip = phase_exhaust_gt(acset, vars)
    expect_identical(vars2flip, exp_out)

    ##cluster_gt
    vars2flip = phase_cluster_gt(acset, vars)
    expect_identical(vars2flip, exp_out)

    ##create acset with ac counts as input
    acset_exp = ac_prephased()
    acset = acset_exp$acset
    exp_out = acset_exp$exp
    vars = rownames(acset[['gt']])
        
    ##wexhaust_gt    
    vars2flip = wphase_exhaust_gt(acset, vars)
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))

    ##wcluster_gt
    vars2flip = wphase_cluster_gt(acset, vars)
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))
})

test_that('all or nothing change when ase have been prephased', {

    ##create acset with ac counts as input
    acset_exp = ac_prephased()
    acset = acset_exp$acset
    exp_out = acset_exp$exp
    vars = rownames(acset[['gt']])

    ##exhaust_ac
    vars2flip = phase_exhaust_ase(acset, vars)
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))

    ##wexhaust_ac
    vars2flip = wphase_exhaust_ase(acset, vars)
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))
    
    ##cluster_ac
    vars2flip = phase_cluster_ase(acset, vars)
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))
    
    ##wcluster_ac
    vars2flip = wphase_cluster_ase(acset, vars)
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))
})

test_that('weights cause a var to be flipped', {

    ##create acset with ac counts as input
    acset_exp = ac_weighedtoflip()
    acset = acset_exp$acset
    exp_out = acset_exp$exp
    vars = rownames(acset[['gt']])    
    
    ##wexhaust_gt
    vars2flip = wphase_exhaust_gt(acset, vars)
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))

    ##wexhaust_ac
    vars2flip = wphase_exhaust_ase(acset, vars)
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))

    ##wcluster_gt
    vars2flip = wphase_cluster_gt(acset, vars)
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))
    
    ##wcluster_ac
    vars2flip = wphase_cluster_ase(acset, vars)
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))
})

test_that('weights cause a var to not be flipped', {

    ##create acset with ac counts as input
    acset_exp = ac_weighedtonotflip()
    acset = acset_exp$acset
    exp_out = acset_exp$exp
    vars = rownames(acset[['gt']])    
    
    ##wexhaust_gt
    vars2flip = wphase_exhaust_gt(acset, vars)
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))

    ##wexhaust_ac
    vars2flip = wphase_exhaust_ase(acset, vars)
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))

    ##wcluster_gt
    vars2flip = wphase_cluster_gt(acset, vars)
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))
    
    ##wcluster_ac
    vars2flip = wphase_cluster_ase(acset, vars)
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))
})
