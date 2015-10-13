

set_aseweights <- function(acset, p = 0.5){

    altcount = acset[['altcount']]
    refcount = acset[['refcount']]

    counts = cbind(as.integer(altcount), as.integer(refcount))

    ##TBD: slow, so either filter, optimize or parallelize
    weights = apply(counts, 1, counts2weight)

    ##vec -> matrix
    weights = matrix(weights, nrow = nrow(altcount))
    rownames(weights) = rownames(altcount)
    colnames(wights) = colnames(altcount)
    
    ##store
    acset[['weights']] = weights

    return(acset)
}

counts2weight <- function(altcount, refcount, p = 0.5){

  weight = 1 - binom.test(altcount, altcount + refcount, p, alternative = 'two.sided')$p.value

  return(weight)
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
    ##TODO: errorcheck that vars and feats has no NA values.
    
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
