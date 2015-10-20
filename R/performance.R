


gtsynth_acset <- function(nvars, nalt_ofvars, ncells, notherhaplo_ofcells, percnoise = 0){
##synthesize genotype matrix and construct acset

    ##haplotypes
    paternal = sample(c(rep(0, nvars - nalt_ofvars), rep(2, nalt_ofvars)))
    maternal = compl_gt(paternal)

    ##create gt matrix
    gt_pat = as.matrix(as.data.frame(rep(list(paternal), ncells - notherhaplo_ofcells)))
    gt_mat = as.matrix(as.data.frame(rep(list(maternal), notherhaplo_ofcells)))
    gt = cbind(gt_pat, gt_mat)
    gt = gt[, sample(1:ncells)]
    vars = 1:nrow(gt)
    colnames(gt) = 1:ncells
    rownames(gt) = vars

    ##featdata
    featdata = as.data.frame(matrix(cbind(rep('jfeat', nvars), 1:nvars), ncol = 2, dimnames = list(vars, c('feat', 'var'))))
    
    ##create acset
    acset = new_acset(featdata, gt = gt)

    return(acset)
}

get_phase_pval <- function(acset, nperm = 100, fixedrowmargin = TRUE){

    ac = !is.null(acset[['refcount']])
    feats = unique(acset[['featdata']][, 'feat'])
    nfeats = length(feats)
    nullscore = matrix(NA, nrow = nfeats, ncol = nperm)
    rownames(nullscore) = feats

    ##get args
    phaseargs = acset[['args']][['phase']]
    input = phaseargs[['input']]
    weigh = phaseargs[['weigh']]
    method = phaseargs[['method']]
    nvars_max = phaseargs[['nvars_max']]
    filterargs = acset[['args']][['filter']]
    nmin_var = filterargs[['nmin_var']]
    min_acount = filterargs[['min_acount']]
    fc = filterargs[['fc']]

    ##*###
    ##Calculate null-dist
    ##*###
    cat('\n')
    progbar = txtProgressBar(min = 0, max = nperm, style = 3)
    for(jperm in 1:nperm){
        
        ##randomize
        if(ac){
            acset_rnd = racset(acset, type = 'ac', fixedrowmargin = fixedrowmargin)
            acset_rnd = filter_nminvar(acset_rnd, nmin_var)
            acset_rnd = call_gt(acset_rnd, min_acount, fc)
        }else{
            acset_rnd = racset(acset, type = 'gt', fixedrowmargin = fixedrowmargin)
        }
        
        ##preproc
        acset_rnd = filter_nminmono(acset_rnd)
        acset_rnd = filter_nminvar(acset_rnd, nmin_var)
        
        ##phase
        acset_rnd = phase(acset_rnd, input, weigh, method, nvars_max, verbosity = 0)

        ##get score
        j.nullscores = acset_rnd[['score']]
        feats = names(j.nullscores)
        nullscore[feats, jperm] = j.nullscores

        setTxtProgressBar(progbar, jperm)
    }


    ##*###
    ##Get p-value
    ##*##
    obsscore = acset[['score']]
    feats = names(obsscore)
    nfeats = length(feats)
    pval = rep(NA, nfeats)
    names(pval) = feats
    
    ##
    for(jfeat in feats){
        pval[jfeat] = (1 + length(which(nullscore[jfeat, ] <= obsscore[jfeat]))) / (nperm + 1)
    }
    
    return(pval)
}

racset <- function(acset, type = 'ac', fixedrowmargin = TRUE){

    if(type == 'ac'){
        acset_rnd = rac(acset, fixedrowmargin)
    }
    if(type == 'gt'){
        acset_rnd = rgt(acset, fixedrowmargin)
    }    
    
    return(acset_rnd)
}

rac <- function(acset, fixedrowmargin = TRUE){

    altcount = acset[['altcount']]
    refcount = acset[['refcount']]
    nrow = nrow(refcount)
    ncol = ncol(refcount)

    if(fixedrowmargin){
        ##permute within row, thus keeping rowmargin fixed
        
        for(jrow in 1:nrow){
            jalt = altcount[jrow, ]
            jref = refcount[jrow, ]
            
            perm.ind = setdiff(1:ncol, which(is.na(jalt) | is.na(jref)))
            perm.sampled.ind = sample(perm.ind)
            jalt[perm.ind] = jalt[perm.sampled.ind]
            jref[perm.ind] = jref[perm.sampled.ind]

            altcount[jrow, ] = jalt
            refcount[jrow, ] = jref
        }
    }else{
        ##Keep table sums of both counts and gts constant. However, row margins are not fixed.

        perm.ind = setdiff(1:(nrow * ncol), which(is.na(altcount) | is.na(refcount)))
        perm.sampled.ind = sample(perm.ind)
        altcount[perm.ind] = altcount[perm.sampled.ind]
        refcount[perm.ind] = refcount[perm.sampled.ind]
    }
    acset_rnd = new_acset(acset[['featdata']], refcount, altcount, acset[['phenodata']])

    return(acset_rnd)
}

rgt <- function(acset, fixedrowmargin = TRUE){
###Some cells may have low read counts, thereby being enriched for NA's. Similarly, they may have high expression, enriching for biallelic calls (gt==1), so not permuting those elements.
    
    gt = acset[['gt']]

    if(fixedrowmargin){
        gt = t(apply(gt, 1, gtperm))
    }else{
        ##Permutation within the whole table, not such that row- or col-margins are fixed.
        gt = gtperm(gt)
    }
    acset[['gt']] = gt

    return(acset)
}

gtperm <- function(gt){
###Some cells may have low read counts, thereby being enriched for NA's. Similarly, they may have high expression, enriching for biallelic calls (gt==1), so not permuting those elements.
    perm.ind = which(gt == 2 | gt == 0)
    perm.sampled.ind = sample(perm.ind)
    gt[perm.ind] = gt[perm.sampled.ind]

    return(gt)
}        

plot_conc <- function(acset, feats = NA, cex = 0.5){

    ##concordance before and after phasing
    conc_pre = acset$gt_conc$conc$feat2ncell
    conc_post = acset$gt_phased_conc$conc$feat2ncell

    ##inconcordance before and after phasing
    inconc_pre = acset$gt_conc$notconc$feat2ncell
    inconc_post = acset$gt_phased_conc$notconc$feat2ncell

    ##order feats
    if(is.logical(feats)){
        conc_post = sort(conc_post, decreasing = TRUE)
        feats = names(conc_post)
    }

    ##ylims
    max_val = max(conc_pre, conc_post, inconc_pre, inconc_post)
    ylim_margin = 1
    ylim_txtmargin = max_val * 0.1
    ylim = c(-max_val - ylim_margin, max_val + ylim_margin)

    ##txt lims
    y_txt = c(max_val - ylim_txtmargin, - max_val + ylim_txtmargin)
    x_txt = c(10, 10)
    
    par(mfrow = c(2, 1), cex = cex)    
    barplot(conc_pre[feats], ylim = ylim, las = 2, names = feats, col = "darkblue", border = "darkblue", main = "Pre-phasing", ylab = "Number of cells")    
    barplot(-inconc_pre[feats], add = TRUE, las = 2, col = "darkgreen", border = "darkgreen")
    text(x = x_txt, y = y_txt, c("Concordant", "Incordordant"), col = c("darkblue", "darkgreen"))
    abline(h = 0)

    barplot(conc_post[feats], ylim = ylim, las = 2, names = feats, col = "darkblue", border = "darkblue", main = "Post-phasing", ylab = "Number of cells")    
    barplot(-inconc_post[feats], add = TRUE, las = 2, col = "darkgreen", border = "darkgreen")
    text(x = x_txt, y = y_txt, c("Concordant", "Incordordant"), col = c("darkblue", "darkgreen"))
    abline(h = 0)
    
}

set_gt_conc <- function(acset){

    featdata = acset[['featdata']]
    acset[['gt_conc']] = calc_gt_conc(acset[['gt']], featdata)
    acset[['gt_phased_conc']] = calc_gt_conc(acset[['gt_phased']], featdata)
    
    return(acset)
}

calc_gt_conc <- function(gt, featdata){

    ##get feats with two vars        
    feat2nvars = table(featdata[, 'feat'])
    feat_pass = names(feat2nvars)[which(feat2nvars == 2)]

    ##subset on vars belonging to feats with two vars
    feat2vars = tapply(featdata[, 'var'], featdata[, 'feat'], unique)
    var_pass = unlist(feat2vars[feat_pass])
    gt = gt[var_pass, ]
    
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

    ##put in list
    concres = list(conc = conc, notconc = notconc)
    
    return(concres)
}

get_conc <- function(feats, samples, inds){
    
    ##feat x sample concordance matrix
    feat2cell = get_concmat(feats, samples, inds)

    ##number of conc cells per feat
    feat2ncells = apply(feat2cell, 1, sum)

    ##put in list
    conc = list(feat2cell = feat2cell, feat2ncells = feat2ncells)

    return(conc)
}

get_concmat <- function(feats, samples, inds){

    nfeats = length(feats)
    nsamples = length(samples)
    feat2cell_conc = matrix(0, nrow = nfeats, ncol = nsamples, dimnames = list(feats, samples))
    feat2cell_conc[inds] = 1
    
    return(feat2cell_conc)
}
