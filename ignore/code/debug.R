

check <- function(){

    ##*###
    ##p-value calc of tcell data with fixed rowsums gives lower p-values than if only fixing the table sum
    ##*###
    ##TBD: an issue specific to {gt, false, cluster}?
    invisible(fakemouse)
    featdata = fakemouse[['featdata']]
    refcount = fakemouse[['refcount']]
    altcount = fakemouse[['altcount']]
    phenodata = fakemouse[['phenodata']]

    acset = new_acset(featdata, refcount, altcount, phenodata)
    lapply(acset, dim)
    acset_rnd = racset(acset)
    lapply(acset, dim)
    nminvar = 2
    acset = filter_nminvar(acset, nminvar)
    lapply(acset, dim)
    min_acount = 3
    fc = 3
    acset = call_gt(acset, min_acount, fc)
    lapply(acset, dim)
    acset = filter_nminmono(acset)
    lapply(acset, dim)
    nminvar = 2
    acset = filter_nminvar(acset, nminvar)
    lapply(acset, dim)

    nperm = 100
    acset = phase(acset, input = 'gt', weigh = FALSE, method = 'cluster')
    pval.row_t = get_phase_pval(acset, nperm, fixedrowmargin = TRUE)
    pval.row_f = get_phase_pval(acset, nperm, fixedrowmargin = FALSE)
    length(which(pval.row_t > 0.99)) #125
    length(which(pval.row_f > 0.99)) #176

    nperm = 10
    acset = phase(acset, input = 'gt', weigh = FALSE, method = 'exhaust')
    pval.row_t = get_phase_pval(acset, nperm, fixedrowmargin = TRUE)
    pval.row_f = get_phase_pval(acset, nperm, fixedrowmargin = FALSE)
    length(which(pval.row_t > 0.99)) #170
    length(which(pval.row_f > 0.99)) #178

    
    ##*###
    ##randomization with fixedrowmargin = T, still reduces the number of monoallelic calls
    ##*###
    setdiff(ac.bak$featdata$feat, acset_rnd$featdata$feat)
    featdata = ac.bak$featdata
    vars = featdata[which(featdata[, 'feat'] == 'ABL2'), 'var']

    gt = acset$gt
    gt.rnd = ac.bak$gt
    j.gt = gt[vars, ]
    j.gt.rnd = gt.rnd[vars, ]
    
    j.gt[which(!is.na(j.gt))]
    j.gt.rnd[which(!is.na(j.gt.rnd))]
    ##second row only contains a single '2'
    
    
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
