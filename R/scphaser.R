
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

racset <- function(acset){
###Randomly swap half of the elements in the count matrices
    
    altcount = acset[['altcount']]
    refcount = acset[['refcount']]

    ##swap half of the elements in the count matrices
    nels = length(altcount)
    swap.ind = sample(1:nels, size = floor(nels / 2))
    
    altcount_swapped = altcount
    refcount_swapped = refcount
    altcount_swapped[swap.ind] = refcount[swap.ind]
    refcount_swapped[swap.ind] = altcount[swap.ind]

    acset_rnd = new_acset(acset[['featdata']], refcount_swapped, altcount_swapped, acset[['phenodata']])
    
    return(acset_rnd)
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

phase <- function(acset, nvars_max = Inf){
    
    ##feat2vars
    featdata = acset[['featdata']]
    feat2var_list = tapply(featdata[, 'var'], featdata[, 'feat'], unique)

    ##record of variants to flip
    varflip = rep(FALSE, nrow(featdata))
    names(varflip) = featdata[, 'var']

    ##loop features
    feats = names(feat2var_list)
    nfeats = length(feats)
    cat('Phasing...\n')
    ##progbar = txtProgressBar(min = 0, max = nfeats, style = 3)
    for(jfeat_it in 1:nfeats){

        jfeat = feats[jfeat_it]
        vars = feat2var_list[[jfeat]]
        nvars = length(vars)        
        
        if(nvars <= nvars_max){
            ##compactness for all possible combinations of flipped vars
            vars2flip = phase_exhaustive(acset, vars)
        }else{
            ##2-class clustering
            vars2flip = phase_cluster(acsest, vars)
        }
        varflip[vars2flip] = TRUE        
        ##setTxtProgressBar(progbar, jfeat_it)
    }
    ##add newline
    cat('\n')

    ##store
    acset[['varflip']] = varflip

    ##add to acset: gt matrix where varflip vars are phased to the complementary gentoype
    acset = set_phased_gt(acset)
    
    return(acset)
}

phase_cluster <- function(acset, vars){

    ##get genotype matrix for vars
    jfeat_gt = acset[['gt']][vars, ]
    jfeat_gt_compl = acset[['gt_compl']][vars, ]
    
    ##cluster
    cluster_ind = cluster::pam(jfeat_gt, k = 2, diss = FALSE, metric = 'euclidean', cluster.only = TRUE)

    ##flip vars belonging to one of the clusters
    vars2flip = vars[which(cluster_ind == 2)]
    gt_jvarflip = jfeat_gt
    gt_jvarflip[vars2flip, ] = jfeat_gt_compl[vars2flip, ]

    ##if the flip caused an increase of compactness then return vars2flip
    flipped_compactness = get_compactness(gt_jvarflip)
    noflip_compactness = get_compactness(jfeat_gt)
    if(flipped_compactness > noflip_compactness){
        vars2flip = vars2flip
    }else{
        vars2flip = character(0)
    }
    
    return(vars2flip)
}

phase_exhaustive <- function(acset, vars){
###calculate a "compactness" score (inversely related to variance) along each cell-axis, and sum up along all cell-axis
###this gives a measure of the compactness of the distribution of all variants in the n-cell space
###our aim is to maximize this compactness        

    ##get genotype matrix for vars
    jfeat_gt = acset[['gt']][vars, ]
    jfeat_gt_compl = acset[['gt_compl']][vars, ]
    
    ##create matrix exhausting all possible combinations of flips
    nvars = length(vars)
    flipcombs = as.matrix(expand.grid(rep(list(c(0, 1)), nvars)))

    ##exclude flipping of all vars since identical results to no change (always last row)
    flipcombs = flipcombs[1:(nrow(flipcombs) - 1), ]

    ##flip a set of variants and calculate compactness
    nflips = nrow(flipcombs)
    flipped_compactness = rep(-1, nflips)
    for(jflip in 1:nflips){

        ##reset gt
        gt_jvarflip = jfeat_gt
        
        ##flip vars2flip to complementary gt
        vars2flip = vars[which(flipcombs[jflip, ] == 1)]
        gt_jvarflip[vars2flip, ] = jfeat_gt_compl[vars2flip, ]
        
        ##compactness of genotype matrix with vars flipped to complementary gt
        flipped_compactness[jflip] = get_compactness(gt_jvarflip)            
    }

    ##get the combination of vars to flip that achieved maximum compactness
    vars2flip = vars[which(flipcombs[which.max(flipped_compactness), ] == 1)]

    return(vars2flip)    
}

phase_singlevarflip <- function(acset){
###TODO: if one wants this to work one needs to update gt_phased within the jvar flip-loop each time the compactness increase.
    
    ##feat2vars
    featdata = acset[['featdata']]
    feat2var_list = tapply(featdata[, 'var'], featdata[, 'feat'], unique)

    ##record of flipped variants
    varflip = rep(FALSE, nrow(featdata))
    names(varflip) = featdata[, 'var']
    
    ##loop features
    feats = names(feat2var_list)
    nfeats = length(feats)
    cat('Phasing...\n')
    progbar = txtProgressBar(min = 0, max = nfeats, style = 3)
    for(jfeat_it in 1:nfeats){

        jfeat = feats[jfeat_it]
        vars = feat2var_list[[jfeat]]
        nvars = length(vars)
        
        ##get genotype matrix
        jfeat_gt = acset[['gt']][vars, ]
        jfeat_gt_compl = acset[['gt_compl']][vars, ]
        gt_phased = jfeat_gt
        
        ##calculate a "compactness" score (inversely related to variance) along each cell-axis, and sum up along all cell-axis
        ##this gives a measure of the compactness of the distribution of all variants in the n-cell space
        ##our aim is to maximize this compactness        
        compactness = get_compactness(jfeat_gt)
        
        ##Flip one variant at a time and check if increasing compactness
        for(jvar in vars){
                            
            jfeat_gt_jvarflip = jfeat_gt
            gt_compl_jvar = jfeat_gt_compl[jvar, ]
            jfeat_gt_jvarflip[jvar, ] = gt_compl_jvar

            ##compactness of genotype matrix with one var flipped to complementary gt
            flipped_compactness = get_compactness(jfeat_gt_jvarflip)

            ##record var if compactness increased
            increase = flipped_compactness > compactness
            if(increase){
                varflip[jvar] = TRUE
                gt_phased[jvar, ] = gt_compl_jvar
            }
        }
        setTxtProgressBar(progbar, jfeat_it)
    }
    ##add newline
    cat('\n')

    ##store
    acset[['varflip']] = varflip

    ##add to acset: gt matrix where varflip vars are phased to the complementary gentoype
    acset = set_phased_gt(acset)
    
    return(acset)
}

get_compactness <- function(gt){
    ##TBD: this also counts elements where only a monoallelic call in a single var. Should be ok, since that number is invariant to a complement shift.
    compactness = sum(abs(colSums(gt == 0, na.rm = TRUE) - colSums(gt == 2, na.rm = TRUE)))

    return(compactness)
}

filter_nminmono <- function(acset, nminmono = 1){
###filter variants having less than nminmono mono-allelic calls of each allele across samples
    
    gt = acset[['gt']]
    pass_vars = which(rowSums(gt == 1, na.rm = TRUE) >= nminmono & rowSums(gt == 2, na.rm = TRUE) >= nminmono)
    acset = subset_rows(acset, pass_vars)

    return(acset)
}

new_acset <- function(featdata, refcount = NA, altcount = NA, phenodata = NA, gt = NA){
    
    ##set phenodata
    if(is.logical(phenodata)){
        if(is.logical(gt)){
            samples = colnames(refcount)
        }else{
            samples = colnames(gt)
        }
        phenodata = as.data.frame(samples, stringsAsFactors = FALSE)
        colnames(phenodata) = 'sample'
        rownames(phenodata) = phenodata[, 'sample']
    }

    ##create data-structure
    if(!is.logical(gt)){
        acset = list(featdata = featdata, phenodata = phenodata, gt = gt)
        acset = set_compl_gt(acset)
    }else{    
        acset = list(featdata = featdata, phenodata = phenodata, refcount = refcount, altcount = altcount)
    }
    
    ##TODO: errorcheck that input featdata and phenodata are dataframes
    ##TODO: check that either refcount and altcount provided and not gt, OR, gt and none of the count matrices
    ##TODO: check that refcount and altcount are matrixes
    ##TODO: check that refcount and altcount have rownames and colnames are not NULL
    ##TODO: errorcheck featdata rownames identical to count matrixes rownames
    ##TODO: errorcheck phenodata rownames identical to count matrixes colnames
    
    return(acset)
}

subset_rows <- function(acset, sel.ind){

    ##data structures with features as rows
    feat_structs = c('altcount', 'refcount', 'featdata')

    if('gt' %in% names(acset)){
        feat_structs = c(feat_structs, 'gt')
    }
    
    if('gt_compl' %in% names(acset)){
        feat_structs = c(feat_structs, 'gt_compl')
    }

    if('gt_phased' %in% names(acset)){
        feat_structs = c(feat_structs, 'gt_phased')
    }    

    ##subset
    acset[feat_structs] = lapply(acset[feat_structs], function(j.obj){j.obj[sel.ind, ]})    

    ##vector data-structures
    if('varflip' %in% names(acset)){
        acset[['varflip']] = acset[['varflip']][sel.ind]
    }        
    
    return(acset)
}

filter_nminvar <- function(acset, nmin_var = 2){
###Filter on minimum number of variants present in feature
    
    ##get n variants per feature
    feat2nvar = table(acset[['featdata']][, 'feat'])

    ##get passed feats
    pass_feat = names(feat2nvar)[which(feat2nvar >= nmin_var)]

    ##subset featdata on passed feats to get passed variants
    featdata = merge(acset[['featdata']], as.matrix(pass_feat), by.x = 'feat', by.y = 1, stringsAsFactors = FALSE)
    rownames(featdata) = featdata[, 'var']
    acset[['featdata']] = featdata
    
    ##subset on passed vars
    pass_var = acset[['featdata']][, 'var']
    acset = subset_rows(acset, pass_var)

    return(acset)
}
##TODO: testthat, identical rownames in all three objects
##TODO: testthat, identical colnames in count mats

call_gt <- function(acset, min_acount = 3, fc = 50){
###Simplistic genotype caller. Genotype callers often rely on DNA-specific assumptions.

    ##indexes of elements with expression
    r.ind = which(acset[['refcount']] >= min_acount)
    a.ind = which(acset[['altcount']] >= min_acount)

    ##Mono-allelic
    ase = acset[['altcount']] / acset[['refcount']]
    a.mono.ind = intersect(a.ind, which(ase >= fc))
    r.mono.ind = intersect(r.ind, which(ase <= 1/fc))

    ##Bi-allelic 
    both.ind = c(a.ind, r.ind)
    both.ind = both.ind[duplicated(both.ind)]
    bi.ind = setdiff(both.ind, c(a.mono.ind, r.mono.ind))

    ##Genotype matrix
    feats = acset[['featdata']][, 'var']
    samples = acset[['phenodata']][, 'sample']
    gt = matrix(data = NA, nrow = length(feats), ncol = length(samples), dimnames = list(feats, samples))
    gt[r.mono.ind] = 0
    gt[bi.ind] = 1
    gt[a.mono.ind] = 2
    
    acset[['gt']] = gt

    ##set complementary genotype matrix
    acset = set_compl_gt(acset)

    return(acset)
}

set_phased_gt <- function(acset){

    gt_phased = acset[['gt']]
    varflip = acset[['varflip']]

    ##phased varflip vars
    gt_phased[varflip, ] = compl_gt(gt_phased[varflip, ])

    ##store
    acset[['gt_phased']] = gt_phased
    
    return(acset)    
}

set_compl_gt <- function(acset){
###Set a complementary genotype matrix where mono-allelic genotypes are swapped

    acset[['gt_compl']] = compl_gt(acset[['gt']])

    return(acset)
}

compl_gt <- function(gt){
###Create complementary genotype matrix where mono-allelic genotypes are swapped

    gt_compl = gt
    gt_compl[which(gt == 0)] = 2
    gt_compl[which(gt == 2)] = 0

    return(gt_compl)
}
