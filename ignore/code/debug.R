

check <- function(){

    ##meth = exhaust, input = ac, weigh = T
    ##slow when nvars: 15
    ##CDC27: nvars: 59
    ##Error in rep.int(rep.int(seq_len(nx), rep.int(rep.fac, nx)), orep) : 
    ##vector is too large
    ##most probably expand.grid with 59 {0,1} binaries.
    
    
    ##*###
    ##Test speed diff btw getting pval by binom.test or pbinom
    ##*###
    altcount = acset[['altcount']]
    refcount = acset[['refcount']]
    counts = cbind(as.integer(altcount), as.integer(refcount))

    ##pbinom
    s = proc.time()
    weights = apply(head(counts, 1e6), 1, counts2weight)    
    diff = proc.time() - s

    ##binom.test
    s = proc.time()
    weights_test = apply(head(counts, 1e6), 1, counts2weight_test)    
    testdiff = proc.time() - s
    ##slower than counts2weight

    
    ##*###
    ##Test: Randomized matrix results differ in concordance after phasing compared to the original
    ##*###

    gt = acset[['gt']]
    gt.rnd = acset_rnd[['gt']]

    which(gt != gt.rnd)
    ##empty.. but, it should differ.

    ##test that randomized allele count mat differs from the input
    altcount = acset[['altcount']]
    altcount.rnd = acset_rnd[['altcount']]
    which(altcount != altcount.rnd)
    ##empty, but should diff.

    refcount = acset[['refcount']]
    refcount.rnd = acset_rnd[['refcount']]
    which(refcount != refcount.rnd)

    
    ##*###
    ##DEBUG: identical concordance before and after phas
    ##*###
    gt = acset[['gt']]
    gt_phased = acset[['gt_phased']]

    featdata = acset[['featdata']]
    jgene = 'APMAP'
    vars = featdata[which(featdata[, 'feat'] == jgene), 'var']
    
    gt = gt[vars, ]
    gt = gt_phased[vars, ]

    var_pass = vars
    
    ##make two matrixes, corresponding to each of the two vars
    nvar = length(var_pass)
    s1_ind = seq(1, nvar-1, by = 2)
    s2_ind = setdiff(1:nvar, s1_ind)
    gt_s1 = gt[s1_ind, ]
    gt_s2 = gt[s2_ind, ]    
    
    ##get elements that have monoallelic calls in both variants
    mono_inds = which((gt_s1 == 0 | gt_s1 == 2) & (gt_s2 == 0 | gt_s2 == 2))

    ##elements where monoallelic genotypes are concordant
    conc_inds = mono_inds[gt_s1[mono_inds] == gt_s2[mono_inds]]
    notconc_inds = setdiff(mono_inds, conc_inds)

    ##feat2cell concordance matrix and feat2ncell counts of concordant gts
    samples = colnames(gt_s1)
    conc = get_conc(feat_pass, samples, conc_inds)
    notconc = get_conc(feat_pass, samples, notconc_inds)
    
}
