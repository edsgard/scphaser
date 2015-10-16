
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

test_that('wphase_exhaustive phases correctly', {

    ##*###
    ##Half of the cells maternal, half paternal
    ##*###

    ##allele count matrixes
    ncells = 10
    ref = c(0, 10, 1, 5, 2, 9)
    alt = c(0, 2, 5, 5, 6, 3)    
    vars2flip_exp = as.integer(c(3, 5))
    
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
    
    ##phase
    vars2flip = wphase_exhaustive(acset, vars)

    ##test
    expect_identical(vars2flip, vars2flip_exp)

    
    ##*###
    ##All cells of one and the same phase
    ##*###
    
    ##allele count matrixes
    ncells = 10
    ref = c(0, 10, 1, 5, 2, 3)
    alt = c(0, 3, 1, 5, 2, 9)    
    vars2flip_exp = as.integer(2)
    
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

    ##phase
    vars2flip = wphase_exhaustive(acset, vars)

    ##test
    expect_identical(vars2flip, vars2flip_exp)
    

    ##*###
    ##Pre-phased (matrix with identical values)
    ##*###

    ##allele count matrixes
    ncells = 10
    ref = c(0, 10, 1, 5, 2, 9)
    alt = c(0, 3, 1, 5, 2, 3)    
    vars2flip_exp = integer(0)
    
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
  
    ##phase
    vars2flip = wphase_exhaustive(acset, vars)

    ##test
    expect_identical(vars2flip, vars2flip_exp)
    
    
    ##*###
    ##Test if weights affect if var should be flipped or not
    ##*###

    ##Weigh the difference heavy, such that a flip should take place
    
    ##allele count matrixes
    refcount = matrix(c(0, 100, 10, 3), nrow = 2)
    altcount = matrix(c(10, 0, 0, 0), nrow = 2)
    vars2flip_exp = as.integer(1)
    
    ncells = ncol(refcount)
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
  
    ##phase
    vars2flip = wphase_exhaustive(acset, vars)

    ##test
    expect_identical(vars2flip, vars2flip_exp)

    
    ##*###
    ##Test if weights affect if var should be flipped or not
    ##*###

    ##Weigh the equal gt heavy, such that no flip should take place
    
    ##allele count matrixes
    refcount = matrix(c(0, 3, 10, 100), nrow = 2)
    altcount = matrix(c(10, 0, 0, 0), nrow = 2)
    vars2flip_exp = as.integer()
    
    ncells = ncol(refcount)
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
  
    ##phase
    vars2flip = wphase_exhaustive(acset, vars)

    ##test
    expect_identical(vars2flip, vars2flip_exp)

    
})

test_that('wphase_cluster_gt phases correctly', {

    ##*###
    ##Half of the cells maternal, half paternal
    ##*###

    ##allele count matrixes
    ncells = 10
    ref = c(0, 10, 1, 5, 2, 9)
    alt = c(0, 2, 5, 5, 6, 3)

    vars2flip_exp = as.integer(c(2, 6))
    
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
    
    ##phase
    vars2flip = wphase_cluster_gt(acset, vars)
    
    ##test
    expect_identical(vars2flip, vars2flip_exp)

    
    ##*###
    ##All cells of one and the same phase
    ##*###
    
    ##allele count matrixes
    ncells = 10
    ref = c(0, 10, 1, 5, 2, 3)
    alt = c(0, 3, 1, 5, 2, 9)    
    vars2flip_exp = as.integer(6)
    
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

    ##phase
    vars2flip = wphase_cluster_gt(acset, vars)

    ##test
    expect_identical(vars2flip, vars2flip_exp)
    

    ##*###
    ##Pre-phased (matrix with identical values)
    ##*###

    ##allele count matrixes
    ncells = 10
    ref = c(0, 10, 1, 5, 2, 9)
    alt = c(0, 3, 1, 5, 2, 3)    
    vars2flip_exp = character(0)
    
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
  
    ##phase
    vars2flip = wphase_cluster_gt(acset, vars)

    ##test
    expect_identical(vars2flip, vars2flip_exp)
    
    
    ##*###
    ##Test if weights affect if var should be flipped or not
    ##*###

    ##Weigh the difference heavy, such that a flip should take place
    
    ##allele count matrixes
    refcount = matrix(c(0, 100, 5, 10, 3, 5), nrow = 3)
    altcount = matrix(c(10, 0, 5, 0, 0, 5), nrow = 3)
    vars2flip_exp = as.integer(2)
    
    ncells = ncol(refcount)
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
  
    ##phase
    vars2flip = wphase_cluster_gt(acset, vars)

    ##test
    expect_identical(vars2flip, vars2flip_exp)

    
    ##*###
    ##Test if weights affect if var should be flipped or not
    ##*###

    ##Weigh the equal gt heavy, such that no flip should take place
    
    ##allele count matrixes
    refcount = matrix(c(0, 3, 5, 10, 100, 5), nrow = 3)
    altcount = matrix(c(10, 0, 5, 0, 0, 5), nrow = 3)
    vars2flip_exp = as.character()
    
    ncells = ncol(refcount)
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
  
    ##phase
    vars2flip = wphase_cluster_gt(acset, vars)

    ##test
    expect_identical(vars2flip, vars2flip_exp)

    
})

test_that('wphase_cluster_ase phases correctly', {

    ##*###
    ##Half of the cells maternal, half paternal
    ##*###

    ##allele count matrixes
    ncells = 10
    ref = c(0, 10, 1, 5, 2, 9)
    alt = c(0, 2, 5, 5, 6, 3)

    vars2flip_exp = as.integer(c(2, 6))
    
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
    
    ##phase
    vars2flip = wphase_cluster_ase(acset, vars)
    
    ##test
    expect_identical(vars2flip, vars2flip_exp)

    
    ##*###
    ##All cells of one and the same phase
    ##*###
    
    ##allele count matrixes
    ncells = 10
    ref = c(0, 10, 1, 5, 2, 3)
    alt = c(0, 3, 1, 5, 2, 9)    
    vars2flip_exp = as.integer(2)
    
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

    ##phase
    vars2flip = wphase_cluster_ase(acset, vars)

    ##test
    expect_identical(vars2flip, vars2flip_exp)
    

    ##*###
    ##Pre-phased (matrix with identical values)
    ##*###

    ##allele count matrixes
    ncells = 10
    ref = c(0, 10, 1, 5, 2, 9)
    alt = c(0, 3, 1, 5, 2, 3)    
    vars2flip_exp = character(0)
    
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
  
    ##phase
    vars2flip = wphase_cluster_ase(acset, vars)

    ##test (none or all flipped)
    expect_true(identical(vars2flip, character(0)) | identical(vars2flip, as.integer(c(2, 6))))
    
    
    ##*###
    ##Test if weights affect if var should be flipped or not
    ##*###

    ##Weigh the difference heavy, such that a flip should take place
    
    ##allele count matrixes
    refcount = matrix(c(0, 100, 5, 10, 3, 5), nrow = 3)
    altcount = matrix(c(10, 0, 5, 0, 0, 5), nrow = 3)
    vars2flip_exp = as.integer(2)
    
    ncells = ncol(refcount)
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
  
    ##phase
    vars2flip = wphase_cluster_ase(acset, vars)

    ##test
    expect_identical(vars2flip, vars2flip_exp)

    
    ##*###
    ##Test if weights affect if var should be flipped or not
    ##*###

    ##Weigh the equal gt heavy, such that no flip should take place
    
    ##allele count matrixes
    refcount = matrix(c(0, 3, 5, 10, 100, 5), nrow = 3)
    altcount = matrix(c(10, 0, 5, 0, 0, 5), nrow = 3)
    vars2flip_exp = as.character()
    
    ncells = ncol(refcount)
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
  
    ##phase
    vars2flip = wphase_cluster_ase(acset, vars)

    ##test
    expect_identical(vars2flip, vars2flip_exp)
    
})
