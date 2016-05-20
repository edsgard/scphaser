
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

    ##expected variants to be flipped: 
    vars2flip_expected = list(as.character(c(2, 5)), as.character(c(1, 3, 4)))
    
    gt = as.matrix(as.data.frame(rep(list(paternal, maternal), ncells / 2)))
    vars = 1:nrow(gt)
    colnames(gt) = 1:ncells
    rownames(gt) = vars

    ##featdata
    nvars = nrow(gt)
    featdata = as.data.frame(matrix(cbind(rep('jfeat', nvars), as.character(1:nvars), rep('dummy', nvars), rep('dummy', nvars)), ncol = 4, dimnames = list(vars, c('feat', 'var', 'ref', 'alt'))), stringsAsFactors = FALSE)
    
    ##create acset
    acset = new_acset(featdata, gt = gt)

    return(list(acset = acset, exp = vars2flip_expected))
}

gt_matpat_twofeat <- function(){

    ##*###
    ##Half of the cells maternal, half paternal
    ##*###

    ##create gt matrix
    ncells = 10
    paternal = c(0, 2, 0, 0, 2)
    maternal = c(2, 0, 2, 2, 0)

    ##expected variants to be flipped: 
    vars2flip_expected = list(jfeat = list(as.character(c(2, 5)), as.character(c(1, 3, 4))), jfeat2 = list(as.character(c(6, 8, 9)), as.character(c(7, 10))))
    
    gt = as.matrix(as.data.frame(rep(list(paternal, maternal), ncells / 2)))
    vars = 1:nrow(gt)
    colnames(gt) = 1:ncells
    rownames(gt) = vars
    
    ##featdata
    nvars = nrow(gt)
    featdata = as.data.frame(matrix(cbind(rep('jfeat', nvars), as.character(1:nvars), rep('dummy', nvars), rep('dummy', nvars)), ncol = 4, dimnames = list(vars, c('feat', 'var', 'ref', 'alt'))), stringsAsFactors = FALSE)

    ##double up
    gt = rbind(gt, gt)
    rownames(gt) = 1:nrow(gt)
    featdata_jfeat2 = featdata
    featdata_jfeat2[, 'feat'] = 'jfeat2'
    featdata = rbind(featdata, featdata_jfeat2)
    featdata[, 'var'] = as.character(1:nrow(featdata))
    rownames(featdata) = featdata[, 'var']
    
    ##create acset
    acset = new_acset(featdata, gt = gt)

    return(list(acset = acset, exp = vars2flip_expected))
}

gt_matpat_hap <- function(){


    ##*###
    ##Half of the cells maternal, half paternal
    ##*###

    ##create gt matrix
    ncells = 10
    paternal = c(0, 2, 0, 0, 2)
    maternal = c(2, 0, 2, 2, 0)

    ##expected variants to be flipped
    vars2flip_expected = list(as.character(c(2, 5)), as.character(c(1, 3, 4)))
    
    gt = as.matrix(as.data.frame(rep(list(paternal, maternal), ncells / 2)))
    vars = 1:nrow(gt)
    colnames(gt) = 1:ncells
    rownames(gt) = vars
    
    ##featdata
    nvars = nrow(gt)
    ref = c('A', 'G', 'G', 'C', 'T')
    alt = c('G', 'C', 'T', 'A', 'C')
    featdata = as.data.frame(matrix(cbind(rep('jfeat', nvars), as.character(1:nvars), ref, alt), ncol = 4, dimnames = list(vars, c('feat', 'var', 'ref', 'alt'))), stringsAsFactors = FALSE)

    ##expected haplotype
    j_exp = vars2flip_expected[[1]]

    vars = rownames(featdata)
    hapA = featdata[, 'ref']
    hapB = featdata[, 'alt']
    names(hapA) = vars
    names(hapB) = vars
    hapA[j_exp] = featdata[j_exp, 'alt']
    hapB[j_exp] = featdata[j_exp, 'ref']
    
    ab = as.data.frame(cbind(hapA, hapB), stringsAsFactors = FALSE)
    ba = ab[, c('hapB', 'hapA')]
    colnames(ba) = colnames(ab)
    hap_exp = list(ab = ab, ba = ba)
    
    ##create acset
    acset = new_acset(featdata, gt = gt)

    res = list(acset = acset, exp = vars2flip_expected, hapexp = hap_exp)
    return(res)
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
    featdata = as.data.frame(matrix(cbind(rep('jfeat', nvars), as.character(1:nvars), rep('dummy', nvars), rep('dummy', nvars)), ncol = 4, dimnames = list(vars, c('feat', 'var', 'ref', 'alt'))), stringsAsFactors = FALSE)

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
    
    ##expected variants to be flipped:
    vars2flip_exp = list(as.character(c(2, 5)), as.character(c(1, 3, 4)))
    
    vars = 1:nrow(gt)
    colnames(gt) = 1:ncells
    rownames(gt) = vars

    ##featdata
    nvars = nrow(gt)    
    featdata = as.data.frame(matrix(cbind(rep('jfeat', nvars), as.character(1:nvars), rep('dummy', nvars), rep('dummy', nvars)), ncol = 4, dimnames = list(vars, c('feat', 'var', 'ref', 'alt'))), stringsAsFactors = FALSE)

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
    featdata = as.data.frame(matrix(cbind(rep('jfeat', nvars), as.character(1:nvars), rep('dummy', nvars), rep('dummy', nvars)), ncol = 4, dimnames = list(vars, c('feat', 'var', 'ref', 'alt'))), stringsAsFactors = FALSE)

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
    
    ##expected variants to be flipped: 
    vars2flip_exp = character(0)
    
    vars = 1:nrow(gt)
    colnames(gt) = 1:ncells
    rownames(gt) = vars

    ##featdata
    nvars = nrow(gt)    
    featdata = as.data.frame(matrix(cbind(rep('jfeat', nvars), as.character(1:nvars), rep('dummy', nvars), rep('dummy', nvars)), ncol = 4, dimnames = list(vars, c('feat', 'var', 'ref', 'alt'))), stringsAsFactors = FALSE)

    
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
    featdata = as.data.frame(matrix(cbind(rep('jfeat', nvars), as.character(1:nvars), rep('dummy', nvars), rep('dummy', nvars)), ncol = 4, dimnames = list(vars, c('feat', 'var', 'ref', 'alt'))), stringsAsFactors = FALSE)

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
    featdata = as.data.frame(matrix(cbind(rep('jfeat', nvars), as.character(1:nvars), rep('dummy', nvars), rep('dummy', nvars)), ncol = 4, dimnames = list(vars, c('feat', 'var', 'ref', 'alt'))), stringsAsFactors = FALSE)

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
    featdata = as.data.frame(matrix(cbind(rep('jfeat', nvars), as.character(1:nvars), rep('dummy', nvars), rep('dummy', nvars)), ncol = 4, dimnames = list(vars, c('feat', 'var', 'ref', 'alt'))), stringsAsFactors = FALSE)

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

ac_onezerocol <- function(){

    refcount = matrix(c(0, 0, 0, 10, 100, 5), nrow = 3)
    altcount = matrix(c(0, 0, 0, 0, 0, 15), nrow = 3)
    vars2flip_exp = list(as.character(3), as.character(c(1, 2)))
    
    ncells = ncol(refcount)
    vars = 1:nrow(refcount)
    samples = 1:ncells
    colnames(refcount) = samples
    rownames(refcount) = vars
    colnames(altcount) = samples
    rownames(altcount) = vars
    
    ##featdata
    nvars = length(vars)
    featdata = as.data.frame(matrix(cbind(rep('jfeat', nvars), as.character(1:nvars), rep('dummy', nvars), rep('dummy', nvars)), ncol = 4, dimnames = list(vars, c('feat', 'var', 'ref', 'alt'))), stringsAsFactors = FALSE)

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

test_that('can phase gt where one cell has all zero counts', {

    ##create gt matrix, providing gt as input
    acset_exp = ac_onezerocol()
    acset = acset_exp$acset
    exp_out = acset_exp$exp
    vars = rownames(acset[['gt']])
    
    ##exhaust_gt
    res = phase_exhaust_gt(acset, vars)
    vars2flip = res$vars2flip
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))

    ##cluster_gt
    res = phase_cluster_gt(acset, vars)
    vars2flip = res$vars2flip
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))

    ##create acset with ac counts as input
    acset_exp = ac_onezerocol()
    acset = acset_exp$acset
    exp_out = acset_exp$exp
    vars = rownames(acset[['gt']])
        
    ##wexhaust_gt
    res = wphase_exhaust_gt(acset, vars)
    vars2flip = res$vars2flip
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))

    ##wcluster_gt
    res = wphase_cluster_gt(acset, vars)
    vars2flip = res$vars2flip
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))            
})

test_that('can phase ase where one cell has all zero counts', {

    ##create gt matrix, providing gt as input
    acset_exp = ac_onezerocol()
    acset = acset_exp$acset
    exp_out = acset_exp$exp
    vars = rownames(acset[['gt']])
    
    ##exhaust_gt
    res = phase_exhaust_ase(acset, vars)
    vars2flip = res$vars2flip
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))

    ##cluster_gt
    res = phase_cluster_ase(acset, vars)
    vars2flip = res$vars2flip
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))

    ##create acset with ac counts as input
    acset_exp = ac_onezerocol()
    acset = acset_exp$acset
    exp_out = acset_exp$exp
    vars = rownames(acset[['gt']])
        
    ##wexhaust_gt
    res = wphase_exhaust_ase(acset, vars)
    vars2flip = res$vars2flip
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))

    ##wcluster_gt
    res = wphase_cluster_ase(acset, vars)
    vars2flip = res$vars2flip
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))            
})

test_that('can phase gt where half cells maternal and half paternal', {

    ##create gt matrix, providing gt as input
    acset_exp = gt_matpat()
    acset = acset_exp$acset
    exp_out = acset_exp$exp
    vars = rownames(acset[['gt']])
    
    ##exhaust_gt
    res = phase_exhaust_gt(acset, vars)
    vars2flip = res$vars2flip
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))

    ##cluster_gt
    res = phase_cluster_gt(acset, vars)
    vars2flip = res$vars2flip
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))

    ##create acset with ac counts as input
    acset_exp = ac_matpat()
    acset = acset_exp$acset
    exp_out = acset_exp$exp
    vars = rownames(acset[['gt']])
        
    ##wexhaust_gt
    res = wphase_exhaust_gt(acset, vars)
    vars2flip = res$vars2flip
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))

    ##wcluster_gt
    res = wphase_cluster_gt(acset, vars)
    vars2flip = res$vars2flip
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))            
})

test_that('haplotype output is correct', {

    ##create gt matrix, providing gt as input
    acset_exp = gt_matpat_hap()
    acset = acset_exp$acset
    exp_out = acset_exp$hapexp
    vars = rownames(acset[['gt']])
    
    ##exhaust_gt
    res = phase(acset, input = 'gt', weigh = FALSE, method = 'exhaust')
    phasedfeat = res$phasedfeat
    hapA = phasedfeat$hapA
    hapB = phasedfeat$hapB
    expect_true((identical(hapA, exp_out[[1]]$hapA) & identical(hapB, exp_out[[1]]$hapB)) | (identical(hapA, exp_out[[2]]$hapA) & identical(hapB, exp_out[[2]]$hapB)))

    ##cluster_gt
    res = phase(acset, input = 'gt', weigh = FALSE, method = 'pam')
    phasedfeat = res$phasedfeat
    hapA = phasedfeat$hapA
    hapB = phasedfeat$hapB
    expect_true((identical(hapA, exp_out[[1]]$hapA) & identical(hapB, exp_out[[1]]$hapB)) | (identical(hapA, exp_out[[2]]$hapA) & identical(hapB, exp_out[[2]]$hapB)))
})

test_that('can phase ase where half cells maternal and half paternal', {

    ##create acset with ac counts as input
    acset_exp = ac_matpat()
    acset = acset_exp$acset
    exp_out = acset_exp$exp
    vars = rownames(acset[['gt']])

    ##exhaust_ac
    res = phase_exhaust_ase(acset, vars)
    vars2flip = res$vars2flip
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))

    ##wexhaust_ac
    res = wphase_exhaust_ase(acset, vars)
    vars2flip = res$vars2flip
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))
    
    ##cluster_ac
    res = phase_cluster_ase(acset, vars)
    vars2flip = res$vars2flip
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))
    
    ##wcluster_ac
    res = wphase_cluster_ase(acset, vars)
    vars2flip = res$vars2flip
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))    
})

test_that('can phase gt where cells have the same phase', {

    ##create gt matrix, providing gt as input
    acset_exp = gt_samephase()
    acset = acset_exp$acset
    exp_out = acset_exp$exp
    vars = rownames(acset[['gt']])
    
    ##exhaust_gt
    res = phase_exhaust_gt(acset, vars)
    vars2flip = res$vars2flip
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))

    ##cluster_gt
    res = phase_cluster_gt(acset, vars)
    vars2flip = res$vars2flip
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))

    ##create acset with ac counts as input
    acset_exp = ac_samephase()
    acset = acset_exp$acset
    exp_out = acset_exp$exp
    vars = rownames(acset[['gt']])
        
    ##wexhaust_gt
    res = wphase_exhaust_gt(acset, vars)
    vars2flip = res$vars2flip
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))

    ##wcluster_gt
    res = wphase_cluster_gt(acset, vars)
    vars2flip = res$vars2flip
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))
})

test_that('can phase ase where cells have the same phase', {

    ##create acset with ac counts as input
    acset_exp = ac_samephase()
    acset = acset_exp$acset
    exp_out = acset_exp$exp
    vars = rownames(acset[['gt']])

    ##exhaust_ac
    res = phase_exhaust_ase(acset, vars)
    vars2flip = res$vars2flip
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))

    ##wexhaust_ac
    res = wphase_exhaust_ase(acset, vars)
    vars2flip = res$vars2flip
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))
    
    ##cluster_ac
    res = phase_cluster_ase(acset, vars)
    vars2flip = res$vars2flip
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))
    
    ##wcluster_ac
    res = wphase_cluster_ase(acset, vars)
    vars2flip = res$vars2flip
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))    
})

test_that('all or nothing change when gt have been prephased', {

    ##create gt matrix, providing gt as input
    acset_exp = gt_prephased()
    acset = acset_exp$acset
    exp_out = acset_exp$exp
    vars = rownames(acset[['gt']])
    
    ##exhaust_gt
    res = phase_exhaust_gt(acset, vars)
    vars2flip = res$vars2flip
    expect_identical(vars2flip, exp_out)

    ##cluster_gt
    res = phase_cluster_gt(acset, vars)
    vars2flip = res$vars2flip
    expect_identical(vars2flip, exp_out)

    ##create acset with ac counts as input
    acset_exp = ac_prephased()
    acset = acset_exp$acset
    exp_out = acset_exp$exp
    vars = rownames(acset[['gt']])
        
    ##wexhaust_gt    
    res = wphase_exhaust_gt(acset, vars)
    vars2flip = res$vars2flip
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))

    ##wcluster_gt
    res = wphase_cluster_gt(acset, vars)
    vars2flip = res$vars2flip
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))
})

test_that('all or nothing change when ase have been prephased', {

    ##create acset with ac counts as input
    acset_exp = ac_prephased()
    acset = acset_exp$acset
    exp_out = acset_exp$exp
    vars = rownames(acset[['gt']])

    ##exhaust_ac
    res = phase_exhaust_ase(acset, vars)
    vars2flip = res$vars2flip
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))

    ##wexhaust_ac
    res = wphase_exhaust_ase(acset, vars)
    vars2flip = res$vars2flip
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))
    
    ##cluster_ac
    res = phase_cluster_ase(acset, vars)
    vars2flip = res$vars2flip
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))
    
    ##wcluster_ac
    res = wphase_cluster_ase(acset, vars)
    vars2flip = res$vars2flip
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))
})

test_that('weights cause a var to be flipped', {

    ##create acset with ac counts as input
    acset_exp = ac_weighedtoflip()
    acset = acset_exp$acset
    exp_out = acset_exp$exp
    vars = rownames(acset[['gt']])    
    
    ##wexhaust_gt
    res = wphase_exhaust_gt(acset, vars)
    vars2flip = res$vars2flip
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))

    ##wexhaust_ac
    res = wphase_exhaust_ase(acset, vars)
    vars2flip = res$vars2flip
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))

    ##wcluster_gt
    res = wphase_cluster_gt(acset, vars)
    vars2flip = res$vars2flip
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))
    
    ##wcluster_ac
    res = wphase_cluster_ase(acset, vars)
    vars2flip = res$vars2flip
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))
})

test_that('weights cause a var to not be flipped', {

    ##create acset with ac counts as input
    acset_exp = ac_weighedtonotflip()
    acset = acset_exp$acset
    exp_out = acset_exp$exp
    vars = rownames(acset[['gt']])    
    
    ##wexhaust_gt
    res = wphase_exhaust_gt(acset, vars)
    vars2flip = res$vars2flip
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))

    ##wexhaust_ac
    res = wphase_exhaust_ase(acset, vars)
    vars2flip = res$vars2flip
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))

    ##wcluster_gt
    res = wphase_cluster_gt(acset, vars)
    vars2flip = res$vars2flip
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))
    
    ##wcluster_ac
    res = wphase_cluster_ase(acset, vars)
    vars2flip = res$vars2flip
    expect_true(identical(vars2flip, exp_out[[1]]) | identical(vars2flip, exp_out[[2]]))
})

test_that('parallelized phasing works', {

    acset_exp = gt_matpat_twofeat()
    acset = acset_exp[['acset']]
    
    ##feat2vars
    featdata = acset[['featdata']]
    feat2var_list = tapply(featdata[, 'var'], featdata[, 'feat'], unique)
    feats = names(feat2var_list)
    nfeats = length(feats)
    
    ##params
    input = 'gt'
    weigh = FALSE
    method = 'exhaust'
    nvars_max = 10
    ##bp_param = BiocParallel::MulticoreParam(workers = 2)
    bp_param = BiocParallel::SerialParam()
    verbosity = 0
    
    ##phase
    res_list = BiocParallel::bplapply(1:nfeats, phase_feat_subset, BPPARAM = bp_param, feat2var_list = feat2var_list, acset = acset, input = input, weigh = weigh, method = method, nvars_max = nvars_max, verbosity = verbosity)
    names(res_list) = feats
    
    ##check results
    feats = names(res_list)
    exp = acset_exp[['exp']]
    for(jfeat in feats){
        jexp = exp[[jfeat]]
        vars2flip = res_list[[jfeat]][['vars2flip']]
        expect_true(identical(vars2flip, jexp[[1]]) | identical(vars2flip, jexp[[2]]))        
    }
    
})
