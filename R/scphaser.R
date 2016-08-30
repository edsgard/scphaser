
#' Phasing of alleles
#'
#' \code{phase} phases alleles using allelic counts or genotype calls from
#' single-cell RNA-seq data.
#'
#' The function phases alleles within each feature. As phasing is not done
#' between features inferred haplotypes should only be used within features.
#'
#' @param acset An acset list created by the \code{\link{new_acset}} function.
#' It must contain elements with either a genotype matrix or two matrixes
#' containing the reference and alternative allele counts, respectively.
#' @param input A character string specifying if allele counts or genotype calls
#' should be used for phasing. Two values are allowed, 'gt' or 'ac'. 'gt'
#' specifies that genotype calls should be used. 'ac' specifies that allele
#' counts should be used.
#' @param weigh A logical specifying if the sample-size, that is, the scale of
#' the allele counts at a variant, should be taken into account. Variants with
#' high counts will be given a greater weight as they are more reliable.
#' @param method A character string specifying the clustering method to be used
#' for the phasing. Two values are allowed, 'exhaust' or 'pam'.
#' @param nvars_max An integer specifying the number of variants within a
#' feature (e.g. gene), above which 'pam' clustering will be used even if the
#' method argument was set to 'exhaust'.
#' @param verbosity An integer specifying the verbosity level. Higher values
#' increase the verbosity.
#' @param bp_param A BiocParallelParam instance, see \code{\link[BiocParallel]{bplapply}}.
#'
#' @return An acset list with the following elements added by the phasing
#' function:
#' 'phasedfeat': Data-frame with six columns. Four first columns, 'feat',
#' 'var', 'ref' and 'alt' are taken from the featdata data-frame of the input
#' acset. The last two columns 'hapA' and 'hapB' contains the haplotype sequence
#' of alleles.
#' 'args': List where each element corresponds to argument values supplied to
#' the phase function.
#' 'varflip': Character vector with names of variants where the alleles were
#' swapped.
#' 'score': Numeric vector with a variability score per feature after phasing.
#' 'gt_phased': Character matrix of genotypes after having swapped alleles
#' according to the inferred phase.
#' 'weights': Numeric matrix with a weight per variant and cell. Only added if
#' arguments set as weigh == TRUE and input == 'gt'.
#'
#' @examples
#' ##create a small artificial genotype matrix
#' ncells = 10
#' paternal = c(0, 2, 0, 0, 2)
#' maternal = c(2, 0, 2, 2, 0)
#' gt = as.matrix(as.data.frame(rep(list(paternal, maternal), ncells / 2)))
#' vars = 1:nrow(gt)
#' colnames(gt) = 1:ncells
#' rownames(gt) = vars
#'
#' ##feature annotation data-frame
#' nvars = nrow(gt)
#' featdata = as.data.frame(matrix(cbind(rep('jfeat', nvars),
#' as.character(1:nvars), rep('dummy', nvars), rep('dummy', nvars)), ncol = 4,
#' dimnames = list(vars, c('feat', 'var', 'ref', 'alt'))), stringsAsFactors =
#' FALSE)
#'
#' ##create acset
#' acset = new_acset(featdata, gt = gt)
#'
#' ##phase
#' acset = phase(acset, input = 'gt', weigh = FALSE, method = 'exhaust',
#' verbosity = 0)
#'
#' ##' The haplotype output is contained in an element of an acset list that was
#' ##added by the phasing function and is named "phasedfeat":
#' head(acset[['phasedfeat']])
#' 
#' @export
phase <- function(acset, input = 'ac', weigh = FALSE, method = 'exhaust', nvars_max = 10, verbosity = -1, bp_param = BiocParallel::SerialParam()){
    
    ##feat2vars
    featdata = acset[['featdata']]
    feat2var_list = tapply(featdata[, 'var'], featdata[, 'feat'], unique)

    ##record of variants to flip
    vars = featdata[, 'var']
    nvars = length(vars)
    varflip = rep(FALSE, nvars)
    names(varflip) = vars

    ##record of score per feature (score = variability between variants after phasing)
    feats = names(feat2var_list)
    nfeats = length(feats)
    score = rep(NA, nfeats)
    names(score) = feats

    ##set weights
    if(weigh & input == 'gt'){
        acset = set_aseweights(acset)
    }
    
    ##loop features
    verbose = R.utils::Verbose(threshold = verbosity)
    R.utils::cat(verbose, '\nPhasing...\n')
    res_list = BiocParallel::bplapply(1:nfeats, phase_feat_subset, BPPARAM = bp_param, feat2var_list = feat2var_list, acset = acset, input = input, weigh = weigh, method = method, nvars_max = nvars_max, verbosity = verbosity)
    names(res_list) = feats
    
    ##store varflip and score
    varflip[unlist(lapply(res_list, '[[', 'vars2flip'))] = TRUE
    score = unlist(lapply(res_list, '[[', 'stat'))
    acset[['varflip']] = names(varflip)[which(varflip)]
    acset[['score']] = score

    ##store haplotype in terms of ref and alt alleles
    acset = set_haplotype(acset)
    
    ##store gt matrix where varflip vars are phased to the complementary gentoype
    acset = set_phased_gt(acset)

    ##store phase arguments
    args = list(input = input, weigh = weigh, method = method, nvars_max = nvars_max)
    acset[['args']][['phase']][c('input', 'weigh', 'method', 'nvars_max')] = args
    
    return(acset)
}

phase_feat_subset <- function(jfeat_it, feat2var_list, acset, input, weigh, method, nvars_max, verbosity = -1){

    ##print progress
    verbose = R.utils::Verbose(threshold = verbosity)
    R.utils::cat(verbose, jfeat_it)
    
    ##subset on vars in jfeat
    vars = feat2var_list[[jfeat_it]]
    jacset = subset_rows(acset, vars)
    
    ##phase single feat
    res = phase_feat(jacset, input, weigh, method, nvars_max)
    
    return(res)
}

phase_feat <- function(acset, input = 'ac', weigh = TRUE, method = 'exhaust', nvars_max = Inf){

    vars = acset[['featdata']][, 'var']
    nvars = length(vars)        
    
    if(method == 'pam'){
        ##can't do 2-class PAM clustering on two observations
        nvars_max = 2
    }
    
    if(!weigh & input == 'gt'){
        if(nvars <= nvars_max){
            ##all possible combinations of flipped vars
            res = phase_exhaust_gt(acset, vars)
        }else{
            ##2-class clustering
            res = phase_cluster_gt(acset, vars)
        }
    }
    
    if(weigh & input == 'gt'){
        if(nvars <= nvars_max){
            ##all possible combinations of flipped vars
            res = wphase_exhaust_gt(acset, vars)
        }else{
            ##2-class clustering
            res = wphase_cluster_gt(acset, vars)
        }
    }

    if(!weigh & input == 'ac'){
        if(nvars <= nvars_max){
            ##all possible combinations of flipped vars
            res = phase_exhaust_ase(acset, vars)
        }else{
            ##2-class clustering
            res = phase_cluster_ase(acset, vars)
        }
    }
    if(weigh & input == 'ac'){
        if(nvars <= nvars_max){
            ##all possible combinations of flipped vars
            res = wphase_exhaust_ase(acset, vars)
        }else{
            ##2-class clustering
            res = wphase_cluster_ase(acset, vars)
        }
    }
    
    return(res)
}

phase_cluster_gt <- function(acset, vars){

    ##get genotype matrix for vars
    jfeat_gt = acset[['gt']][vars, ]
    jfeat_gt_compl = acset[['gt_compl']][vars, ]

    ##the distance need to be increased if NAs. Treat it the same as biallelic.
    jfeat_gt[which(is.na(jfeat_gt))] = 1

    ##cluster
    cluster_ind = cluster::pam(jfeat_gt, k = 2, diss = FALSE, metric = 'euclidean', cluster.only = TRUE)

    ##flip vars belonging to one of the clusters
    vars2flip = vars[which(cluster_ind == 2)]
    gt_jvarflip = jfeat_gt
    gt_jvarflip[vars2flip, ] = jfeat_gt_compl[vars2flip, ]

    ##if the flip caused a decrease in variability then return vars2flip
    flipped_variability = get_var_gt(gt_jvarflip)
    noflip_variability = get_var_gt(jfeat_gt)
    if(flipped_variability < noflip_variability){
        vars2flip = vars2flip
        min_var = flipped_variability
    }else{
        vars2flip = character(0)
        min_var = noflip_variability
    }
    
    return(list(vars2flip = vars2flip, stat = min_var))
}

wphase_cluster_gt <- function(acset, vars){
    
    ##get data-structures
    acset = subset_rows(acset, vars)
    jfeat_gt = acset[['gt']]
    jfeat_gt_compl = acset[['gt_compl']]
    altcount = acset[['altcount']]
    refcount = acset[['refcount']]

    ##the distance need to be increased if NAs. Treat it the same as biallelic.
    jfeat_gt[which(is.na(jfeat_gt))] = 1
    
    ##calculate weighted pairwise distance matrix between variants
    ##let weights for variant i, j in cell k correspond to w_ijk = sum(counts(2x2 table_ij)) / sum_k(2x2_ij tables)
    nvars = length(vars)
    distmat = matrix(NA, nrow = nvars, ncol = nvars)
    for(jvar_it in 1:(nvars-1)){
        for(kvar_it in (jvar_it + 1):nvars){

            ##weights
            jvars = c(jvar_it, kvar_it)
            jvars_count = apply(altcount[jvars, ] + refcount[jvars, ], 2, sum, na.rm = TRUE)
            weights = jvars_count / sum(jvars_count)

            ##weighted Manhattan distance            
            distmat[kvar_it, jvar_it] = sum(weights * abs(jfeat_gt[jvar_it, ] - jfeat_gt[kvar_it, ]))                
        }
    }
    
    ##cluster
    cluster_ind = cluster::pam(stats::as.dist(distmat), k = 2, diss = TRUE, cluster.only = TRUE)

    ##flip vars belonging to one of the clusters
    vars2flip = vars[which(cluster_ind == 2)]

    ##if the flip caused a decrease of the variability then return vars2flip
    gt_jvarflip = jfeat_gt
    gt_jvarflip[vars2flip, ] = jfeat_gt_compl[vars2flip, ]

    ##variability
    weights = acset[['weights']]
    flipped_variability = get_wvar_gt(gt_jvarflip, weights)
    noflip_variability = get_wvar_gt(jfeat_gt, weights)
    if(flipped_variability < noflip_variability){
        vars2flip = vars2flip
        min_var = flipped_variability
    }else{
        vars2flip = character(0)
        min_var = noflip_variability
    }
    
    return(list(vars2flip = vars2flip, stat = min_var))
}

phase_cluster_ase <- function(acset, vars){
    
    ##get data-structures
    acset = subset_rows(acset, vars)
    altcount = acset[['altcount']]
    refcount = acset[['refcount']]

    ##ase
    totcount = altcount + refcount
    ase = (altcount) / totcount

    ##the distance need to be increased if NAs. Treat it the same as biallelic.
    ##handles 0/0 division where both alt and ref counts are 0.
    ase[is.na(ase)] = 0.5        
    
    ##calculate distance matrix between pairs of variants
    var_dist = stats::dist(ase, method = 'euclidean')
    
    ##cluster
    cluster_ind = cluster::pam(var_dist, k = 2, diss = TRUE, cluster.only = TRUE)

    ##flip vars belonging to one of the clusters
    vars2flip = vars[which(cluster_ind == 2)]

    ##if the flip caused a decrease in variability then return vars2flip
    ##flip 
    ase_jvarflip = ase
    ase_compl = 1 - ase
    ase_jvarflip[vars2flip, ] = ase_compl[vars2flip, ]
    
    ##variability
    ##flipped_cv = get_cv_ase(ase_jvarflip)
    ##noflip_cv = get_cv_ase(ase)
    flipped_variability = get_var_ase(ase_jvarflip)
    noflip_variability = get_var_ase(ase)

    if(flipped_variability < noflip_variability){
        vars2flip = vars2flip
        min_var = flipped_variability
    }else{
        vars2flip = character(0)
        min_var = noflip_variability
    }
    
    return(list(vars2flip = vars2flip, stat = min_var))
}

wphase_cluster_ase <- function(acset, vars){
    
    ##get data-structures
    acset = subset_rows(acset, vars)
    altcount = acset[['altcount']]
    refcount = acset[['refcount']]

    ##ase
    totcount = altcount + refcount
    ase = (altcount) / totcount

    ##the distance need to be increased if NAs. Treat it the same as biallelic.
    ##handles 0/0 division where both alt and ref counts are 0.
    ase[is.na(ase)] = 0.5        
    
    ##calculate weighted pairwise distance matrix between variants
    ##let weights for variant i, j in cell k correspond to w_ijk = sum(counts(2x2 table_ij)) / sum_k(2x2_ij tables)
    nvars = length(vars)
    distmat = matrix(NA, nrow = nvars, ncol = nvars)
    for(jvar_it in 1:(nvars-1)){
        for(kvar_it in (jvar_it + 1):nvars){

            ##weights
            jvars = c(jvar_it, kvar_it)
            jvars_count = apply(altcount[jvars, ] + refcount[jvars, ], 2, sum, na.rm = TRUE)
            weights = jvars_count / sum(jvars_count)

            ##weighted distance
            distmat[kvar_it, jvar_it] = sum(weights * abs(ase[jvar_it, ] - ase[kvar_it, ]))
            
            ##c = boot::corr(t(jfeat_ase[jvars, ]), w = weights)            
            ##d = ((-1 * c) / 2) + 0.5            
        }
    }
    
    ##cluster
    cluster_ind = cluster::pam(stats::as.dist(distmat), k = 2, diss = TRUE, cluster.only = TRUE)

    ##flip vars belonging to one of the clusters
    vars2flip = vars[which(cluster_ind == 2)]

    ##if the flip caused a decrease in variability then return vars2flip
    ##flip 
    ase_jvarflip = ase
    ase_compl = 1 - ase
    ase_jvarflip[vars2flip, ] = ase_compl[vars2flip, ]
    
    ##variability
    weights = totcount
    ##flipped_cv = get_wcv_ase(ase_jvarflip, weights)
    ##noflip_cv = get_wcv_ase(ase, weights)
    flipped_variability = get_wvar_ase(ase_jvarflip, weights)
    noflip_variability = get_wvar_ase(ase, weights)    
    if(flipped_variability < noflip_variability){
        vars2flip = vars2flip
        min_var = flipped_variability
    }else{
        vars2flip = character(0)
        min_var = noflip_variability
    }

    return(list(vars2flip = vars2flip, stat = min_var))
}

phase_exhaust_gt <- function(acset, vars){
###calculate a "compactness" score along each cell-axis, and sum up along all cell-axis to get a total compactness score of the distribution of variants.
###1 / compactness gives a measure of the variability of the distribution of all variants in the n-cell space
###our aim is to minimize the variablity

    ##get genotype matrix for vars
    jfeat_gt = acset[['gt']][vars, ]
    jfeat_gt_compl = acset[['gt_compl']][vars, ]
    
    ##create matrix exhausting all possible combinations of flips
    nvars = length(vars)
    flipcombs = as.matrix(expand.grid(rep(list(c(0, 1)), nvars)))

    ##exclude flipping of all vars since identical results to no change (always last row)
    flipcombs = flipcombs[1:(nrow(flipcombs) - 1), ]

    ##flip a set of variants and calculate variability
    nflips = nrow(flipcombs)
    flipped_variability = rep(Inf, nflips)
    for(jflip in 1:nflips){

        ##reset gt
        gt_jvarflip = jfeat_gt
        
        ##flip vars2flip to complementary gt
        vars2flip = vars[which(flipcombs[jflip, ] == 1)]
        gt_jvarflip[vars2flip, ] = jfeat_gt_compl[vars2flip, ]
        
        ##variability of genotype matrix with vars flipped to complementary gt
        flipped_variability[jflip] = get_var_gt(gt_jvarflip)            
    }

    ##get the combination of vars to flip that achieved minimum variability
    min_ind = which.min(flipped_variability)
    vars2flip = vars[which(flipcombs[min_ind, ] == 1)]
    min_variability = flipped_variability[min_ind]
    
    return(list(vars2flip = vars2flip, stat = min_variability))
}

wphase_exhaust_gt <- function(acset, vars){

    ##get genotype and weight matrix for vars
    jfeat_gt = acset[['gt']][vars, ]
    jfeat_gt_compl = acset[['gt_compl']][vars, ]
    weights = acset[['weights']][vars, ]
    
    ##create matrix exhausting all possible combinations of flips
    nvars = length(vars)
    flipcombs = as.matrix(expand.grid(rep(list(c(0, 1)), nvars)))

    ##exclude flipping of all vars since identical results to no change (always last row)
    flipcombs = flipcombs[1:(nrow(flipcombs) - 1), ]

    ##flip a set of variants and calculate variability
    nflips = nrow(flipcombs)
    flipped_variability = rep(Inf, nflips)
    for(jflip in 1:nflips){

        ##reset gt
        gt_jvarflip = jfeat_gt
        
        ##flip vars2flip to complementary gt
        vars2flip = vars[which(flipcombs[jflip, ] == 1)]
        gt_jvarflip[vars2flip, ] = jfeat_gt_compl[vars2flip, ]
        
        ##weighted variability of genotype matrix with vars flipped to complementary gt
        flipped_variability[jflip] = get_wvar_gt(gt_jvarflip, weights)            
    }

    ##get the combination of vars to flip that achieved maximum variability
    min_ind = which.min(flipped_variability)
    vars2flip = vars[which(flipcombs[min_ind, ] == 1)]
    min_variability = flipped_variability[min_ind]

    return(list(vars2flip = vars2flip, stat = min_variability))
}

phase_exhaust_ase <- function(acset, vars){

    ##get data-structures
    acset = subset_rows(acset, vars)
    altcount = acset[['altcount']]
    refcount = acset[['refcount']]

    ##ase
    totcount = altcount + refcount
    ase = altcount / totcount

    ##the distance need to be increased if NAs. Treat it the same as biallelic.
    ##handles 0/0 division where both alt and ref counts are 0.
    ase[is.na(ase)] = 0.5

    ##complement
    ase_compl = 1 - ase
    
    ##create matrix exhausting all possible combinations of flips
    nvars = length(vars)
    flipcombs = as.matrix(expand.grid(rep(list(c(0, 1)), nvars)))

    ##exclude flipping of all vars since identical results to no change (always last row)
    flipcombs = flipcombs[1:(nrow(flipcombs) - 1), ]

    ##flip a set of variants and calculate variability
    nflips = nrow(flipcombs)
    flipped_variability = rep(Inf, nflips)
    for(jflip in 1:nflips){

        ##reset ase matrix
        ase_jvarflip = ase
        
        ##flip vars2flip to complementary ase
        vars2flip = vars[which(flipcombs[jflip, ] == 1)]
        ase_jvarflip[vars2flip, ] = ase_compl[vars2flip, ]
        
        ##cv of ase matrix with vars flipped to complementary ase
        ##flipped_cv[jflip] = get_cv_ase(ase_jvarflip)
        flipped_variability[jflip] = get_var_ase(ase_jvarflip)
    }

    ##get the combination of vars to flip that achieved minimum cv
    min_ind = which.min(flipped_variability)
    vars2flip = vars[which(flipcombs[min_ind, ] == 1)]    
    min_var = flipped_variability[min_ind]
   
    return(list(vars2flip = vars2flip, stat = min_var))
}

wphase_exhaust_ase <- function(acset, vars){

    ##get data-structures
    acset = subset_rows(acset, vars)
    altcount = acset[['altcount']]
    refcount = acset[['refcount']]

    ##ase
    totcount = altcount + refcount
    ase = altcount / totcount

    ##the distance need to be increased if NAs. Treat it the same as biallelic.
    ##handles 0/0 division where both alt and ref counts are 0.
    ase[is.na(ase)] = 0.5        

    ##complement
    ase_compl = 1 - ase
    
    ##create matrix exhausting all possible combinations of flips
    nvars = length(vars)
    flipcombs = as.matrix(expand.grid(rep(list(c(0, 1)), nvars)))

    ##exclude flipping of all vars since identical results to no change (always last row)
    flipcombs = flipcombs[1:(nrow(flipcombs) - 1), ]

    ##flip a set of variants and calculate variability
    nflips = nrow(flipcombs)
    flipped_variability = rep(Inf, nflips)
    for(jflip in 1:nflips){

        ##reset ase matrix
        ase_jvarflip = ase
        
        ##flip vars2flip to complementary ase
        vars2flip = vars[which(flipcombs[jflip, ] == 1)]
        ase_jvarflip[vars2flip, ] = ase_compl[vars2flip, ]
        
        ##cv of ase matrix with vars flipped to complementary ase
        ##flipped_cv[jflip] = get_wcv_ase(ase_jvarflip, weights = totcount)
        flipped_variability[jflip] = get_wvar_ase(ase_jvarflip, weights = totcount)
    }

    ##get the combination of vars to flip that achieved minimum cv
    min_ind = which.min(flipped_variability)
    vars2flip = vars[which(flipcombs[min_ind, ] == 1)]    
    min_var = flipped_variability[min_ind]
   
    return(list(vars2flip = vars2flip, stat = min_var))
}

get_var_gt_mat <- function(acset){
    
    ##feat2vars
    featdata = acset[['featdata']]
    feat2var_list = tapply(featdata[, 'var'], featdata[, 'feat'], unique)

    ##loop feats    
    feats = names(feat2var_list)
    score = rep(NA, length(feats))
    names(score) = feats
    gt = acset[['gt']]
    for(jfeat in feats){
        vars = featdata[which(featdata[, 'feat'] == jfeat), 'var']
        jgt = gt[vars, ]
        score[jfeat] = get_var_gt(jgt)
    }

    return(score)
}

get_var_gt <- function(gt){
    ##TBD: this also counts elements where only a monoallelic call in a single var. Should be ok, since that number is invariant to a complement shift.
    compactness = sum(abs(colSums(gt == 0, na.rm = TRUE) - colSums(gt == 2, na.rm = TRUE)))

    v = 1 / compactness
    
    return(v)
}

get_wvar_gt <- function(gt, weights){
###weighted compactness

    nvars = nrow(gt)
    nsamples = ncol(gt)
    
    ref_mat = matrix(0, nrow = nvars, ncol = nsamples)    
    ref_ind = which(gt == 0)
    ref_mat[ref_ind] = weights[ref_ind]
    
    alt_mat = matrix(0, nrow = nvars, ncol = nsamples)
    alt_ind = which(gt == 2)
    alt_mat[alt_ind] = weights[alt_ind]
    
    compactness = sum(abs(colSums(ref_mat) - colSums(alt_mat)))

    v = 1 / compactness
    
    return(v)
}

get_var_ase <- function(ase){

    ##sum the vars from every cell (assume cells are independent observations)
    vartot = sum(apply(ase, 2, stats::var), na.rm = TRUE)

    return(vartot)
}

get_wvar_ase <- function(ase, weights){

    ##sum the weighted vars from every cell (assume cells are independent observations)
    ncells = ncol(ase)
    vars = rep(NA, ncells)
    for(jcell in 1:ncells){
        jase = ase[, jcell]
        jw = weights[, jcell]
        vars[jcell] = Hmisc::wtd.var(jase, jw)
    }
    vartot = sum(vars, na.rm = TRUE)

    return(vartot)
}

set_aseweights <- function(acset, p = 0.5){

    if(is.null(acset[['altcount']]) || is.null(acset[['refcount']])){
        stop('altcount or refcount is null')
    }
    
    altcount = acset[['altcount']]
    refcount = acset[['refcount']]
    
    ##matrix -> vec for fast calculation
    counts = cbind(as.integer(altcount), as.integer(refcount))
    counts = as.data.frame(counts)

    ##set weight to 0 if altcount == refcount
    weights = rep(0, nrow(counts))
    diff.ind = which(counts[, 1] != counts[, 2])
    
    ##calculate weight as 1 - two.sided p-value
    weights[diff.ind] = apply(counts[diff.ind, ], 1, counts2weight)
    
    ##vec -> matrix
    weights = matrix(weights, nrow = nrow(altcount))
    rownames(weights) = rownames(altcount)
    colnames(weights) = colnames(altcount)
    
    ##store
    acset[['weights']] = weights

    return(acset)
}

counts2weight <- function(counts, p = 0.5){
###1 - (two-sided p-value)
###Using pbinom is faster than binom.test (weight = 1 - binom.test(counts[1], sum(counts), p, alternative = 'two.sided')$p.value)
    
    if(counts[1] == counts[2]){
        weight = 0
    }else{
        weight = 1 - 2 * stats::pbinom(min(counts), sum(counts), p)
    }
    return(weight)
}

#' Construct an allele count set
#'
#' \code{new_acset} constructs a list which is the data-structure used by the
#' scphaser phasing functions.
#'
#' The function performs a number of error-checks to ensure that the constructed
#' acset-list elements satisfy the data-format used by the phasing functions.
#'
#' @param featdata Data-frame with four required columns and arbitrary
#' additional columns. The purpose of the featdata data-frame is to map variants
#' to features and specify the two alleles of each variant. The four required
#' columns must be named 'feat', 'var', 'ref' and 'alt'. The rownames must also
#' be set to be identical to the var column. 'feat' is a character vector
#' specifying feature names, such as gene names. 'var' is a character vector
#' specifying variant names, such as dbSNP rs id. 'ref' and 'alt' are character
#' vectors specifying the alleles of each variant, such as the reference and
#' alternative allele.
#' @param refcount Matrix with allelic counts for the reference allele with
#' variants as rows and cells as columns. The rownames have to match the values
#' in the 'var' column in the featdata, and the colnames the values in the
#' phenodata 'sample' column. If refcount is provided altcount must also be
#' provided. If the gt argument is provided then refcount and altcount are not
#' required arguments.
#' @param altcount Matrix with allelic counts for the alternative allele, where
#' the row- and col-names must match those of refcount.
#' @param gt Matrix with integer values representing transcribed genotype calls.
#' 0: reference allele most highly expressed, 1: bi-allelic expression with
#' similar degree of expression from the two alleles, 2: alternative allele most
#' highly expressed. NA's are allowed and can be used to indicate entries where
#' no call could be made. The rownames have to match the values in the 'var'
#' column in the featdata, and the colnames the values in the phenodata 'sample'
#' column. If refcount and altcount are provided then the gt argument does not
#' need to be provided.
#' @param phenodata Data-frame which annotates the cells. It must contain a
#' column named 'sample'. If the phenodata argument is not provided it is
#' created by the \code{link{new_acset}} function with the sample column set to
#' be identical to the column names of refcount or gt.
#' 
#' @return acset An acset list which contains elements required to apply the
#' phasing functions.
#'
#' @examples
#' ##create a small artificial genotype matrix
#' ncells = 10
#' paternal = c(0, 2, 0, 0, 2)
#' maternal = c(2, 0, 2, 2, 0)
#' gt = as.matrix(as.data.frame(rep(list(paternal, maternal), ncells / 2)))
#' vars = 1:nrow(gt)
#' colnames(gt) = 1:ncells
#' rownames(gt) = vars
#'
#' ##create a feature annotation data-frame
#' nvars = nrow(gt)
#' featdata = as.data.frame(matrix(cbind(rep('jfeat', nvars),
#' as.character(1:nvars), rep('dummy', nvars), rep('dummy', nvars)), ncol = 4,
#' dimnames = list(vars, c('feat', 'var', 'ref', 'alt'))), stringsAsFactors =
#' FALSE)
#'
#' ##create acset
#' acset = new_acset(featdata, gt = gt)
#' 
#' @export
new_acset <- function(featdata, refcount = NA, altcount = NA, phenodata = NA, gt = NA){

    ##Error-checks
    if(is.logical(gt) && (is.logical(refcount) || is.logical(altcount))){
        stop('Neither gt nor refcount and altcount has been provided')
    }
    if(!is.logical(refcount) && is.logical(altcount)){
        stop('refcount provided but not altcount')
    }
    if(is.logical(refcount) && !is.logical(altcount)){
        stop('altcount provided but not refcount')
    }
    
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
        if(is.logical(refcount)){
            acset = list(featdata = featdata, phenodata = phenodata, gt = gt)
        }else{
            acset = list(featdata = featdata, phenodata = phenodata, gt = gt, refcount = refcount, altcount = altcount)
        }
        acset = set_compl_gt(acset)
    }else{    
        acset = list(featdata = featdata, phenodata = phenodata, refcount = refcount, altcount = altcount)
    }

    ##*###
    ##Error-checks
    ##*###
    if(!('var' %in% colnames(featdata))){
        stop('featdata require a column named "var"')
    }
    if(!('feat' %in% colnames(featdata))){
        stop('featdata require a column named "feat"')
    }
    if(!('ref' %in% colnames(featdata))){
        stop('featdata require a column named "ref"')
    }
    if(!('alt' %in% colnames(featdata))){
        stop('featdata require a column named "alt"')
    }
    if(!('sample' %in% colnames(phenodata))){
        stop('phenodata require a column named "sample"')
    }
    if(!(is.character(featdata[, 'var']))){
        stop('The column "var" in featdata need to be of class "character"')
    }
    if(!(is.character(featdata[, 'feat']))){
        stop('The column "feat" in featdata need to be of class "character"')
    }
    if(!(is.character(featdata[, 'ref']))){
        stop('The column "ref" in featdata need to be of class "character"')
    }
    if(!(is.character(featdata[, 'alt']))){
        stop('The column "alt" in featdata need to be of class "character"')
    }
    if(!(is.data.frame(featdata))){
        stop('featdata need to be a data-frame')
    }
    if(!(is.data.frame(phenodata))){
        stop('phenodata need to be a data-frame')
    }
    if(length(which(is.na(featdata[, 'var']))) > 0){
        stop('NA values in "var" column of featdata')
    }
    if(length(which(is.na(featdata[, 'feat']))) > 0){
        stop('NA values in "feat" column of featdata')
    }
                  
    if(!is.logical(refcount)){
        if(!(is.matrix(refcount))){
            stop('refcount need to be a matrix')
        }
        if(!(is.matrix(altcount))){
            stop('altcount need to be a matrix')
        }
        if(is.null(rownames(refcount))){
            stop('rownames of refcount is null')
        }
        if(is.null(rownames(altcount))){
            stop('rownames of altcount is null')
        }
        if(is.null(colnames(refcount))){
            stop('colnames of refcount is null')
        }
        if(is.null(colnames(altcount))){
            stop('colnames of altcount is null')
        }
        if(!(identical(rownames(featdata), rownames(altcount)) && identical(rownames(featdata), rownames(refcount)))){
            stop('rownames of featdata and refcount and altcount need to be identical')
        }
        if(!(identical(rownames(phenodata), colnames(altcount)) && identical(rownames(phenodata), colnames(refcount)))){
            stop('rownames of phenodata need to be identical to the colnames of refcount and altcount')
        }
    }
    if(!is.logical(gt)){
        if(length(setdiff(gt, c(NA, 0, 1, 2)))){
            stop('gt contains values that are not in the allowed set {0, 1, 2, NA}')
        }
        if(!(is.matrix(gt))){
            stop('gt need to be a matrix')
        }
        if(is.null(rownames(gt))){
            stop('rownames of gt is null')
        }
        if(is.null(colnames(gt))){
            stop('colnames of gt is null')
        }
        if(!(identical(rownames(featdata), rownames(gt)))){
            stop('rownames of featdata and gt need to be identical')
        }
        if(!(identical(rownames(phenodata), colnames(gt)))){
            stop('rownames of phenodata need to be identical to the colnames of gt')
        }        
    }
    
    return(acset)
}

#' Subset variants
#'
#' \code{subset_rows} subsets the variants of an acset
#'
#' @param acset An acset list created by the \code{\link{new_acset}} function.
#' @param sel.ind A vector with either integer row indices or a character
#' vector to subset the rownames.
#' 
#' @return acset An acset list where every element containing variant-related
#' information has been subsetted
#'
#' @examples
#' ##create a small artificial genotype matrix
#' ncells = 10
#' paternal = c(0, 2, 0, 0, 2)
#' maternal = c(2, 0, 2, 2, 0)
#' gt = as.matrix(as.data.frame(rep(list(paternal, maternal), ncells / 2)))
#' vars = 1:nrow(gt)
#' colnames(gt) = 1:ncells
#' rownames(gt) = vars
#'
#' ##create a feature annotation data-frame
#' nvars = nrow(gt)
#' featdata = as.data.frame(matrix(cbind(rep('jfeat', nvars),
#' as.character(1:nvars), rep('dummy', nvars), rep('dummy', nvars)), ncol = 4,
#' dimnames = list(vars, c('feat', 'var', 'ref', 'alt'))), stringsAsFactors =
#' FALSE)
#'
#' ##create acset
#' acset = new_acset(featdata, gt = gt)
#'
#' ##subset variants
#' sel_ind = c(1, 3)
#' acset = subset_rows(acset, sel_ind)
#' 
#' @export
subset_rows <- function(acset, sel.ind){

    ##data structures with features as rows
    feat_structs = c('featdata')

    if('refcount' %in% names(acset)){
        feat_structs = c(feat_structs, 'refcount')
    }

    if('altcount' %in% names(acset)){
        feat_structs = c(feat_structs, 'altcount')
    }

    if('gt' %in% names(acset)){
        feat_structs = c(feat_structs, 'gt')
    }
    
    if('gt_compl' %in% names(acset)){
        feat_structs = c(feat_structs, 'gt_compl')
    }

    if('gt_phased' %in% names(acset)){
        feat_structs = c(feat_structs, 'gt_phased')
    }    

    if('weights' %in% names(acset)){
        feat_structs = c(feat_structs, 'weights')
    }    

    if('phasedfeat' %in% names(acset)){
        feat_structs = c(feat_structs, 'phasedfeat')
    }        
    
    ##subset
    acset[feat_structs] = lapply(acset[feat_structs], function(j.obj){j.obj[sel.ind, ]})    

    ##vector data-structures
    if('varflip' %in% names(acset)){
        acset[['varflip']] = intersect(acset[['varflip']], acset[['featdata']][, 'var'])
    }        
    
    return(acset)
}

#' Subset cells
#'
#' \code{subset_cols} subsets the cells of an acset
#'
#' @param acset An acset list created by the \code{\link{new_acset}} function.
#' @param sel.ind A vector with either integer column indices or a character
#' vector to subset the colnames.
#' 
#' @return acset An acset list where every element containing cell-related
#' information has been subsetted
#'
#' @examples
#' ##create a small artificial genotype matrix
#' ncells = 10
#' paternal = c(0, 2, 0, 0, 2)
#' maternal = c(2, 0, 2, 2, 0)
#' gt = as.matrix(as.data.frame(rep(list(paternal, maternal), ncells / 2)))
#' vars = 1:nrow(gt)
#' colnames(gt) = 1:ncells
#' rownames(gt) = vars
#'
#' ##create a feature annotation data-frame
#' nvars = nrow(gt)
#' featdata = as.data.frame(matrix(cbind(rep('jfeat', nvars),
#' as.character(1:nvars), rep('dummy', nvars), rep('dummy', nvars)), ncol = 4,
#' dimnames = list(vars, c('feat', 'var', 'ref', 'alt'))), stringsAsFactors =
#' FALSE)
#'
#' ##create acset
#' acset = new_acset(featdata, gt = gt)
#'
#' ##subset variants
#' sel_ind = c(1, 3)
#' acset = subset_cols(acset, sel_ind)
#' 
#' @export
subset_cols <- function(acset, sel.ind){

    ##data structures with features as rows
    acset[['phenodata']] = acset[['phenodata']][sel.ind, , drop = FALSE]

    pheno_structs = NULL
    if('refcount' %in% names(acset)){
        pheno_structs = c(pheno_structs, 'refcount')
    }

    if('altcount' %in% names(acset)){
        pheno_structs = c(pheno_structs, 'altcount')
    }

    if('gt' %in% names(acset)){
        pheno_structs = c(pheno_structs, 'gt')
    }
    
    if('gt_compl' %in% names(acset)){
        pheno_structs = c(pheno_structs, 'gt_compl')
    }

    if('gt_phased' %in% names(acset)){
        pheno_structs = c(pheno_structs, 'gt_phased')
    }    

    ##subset
    acset[pheno_structs] = lapply(acset[pheno_structs], function(j.obj){j.obj[, sel.ind]})
    
    return(acset)
}

#' Subset features
#'
#' \code{subset_feat} subsets the features of an acset
#'
#' @param acset An acset list created by the \code{\link{new_acset}} function.
#' @param feat A character vector with features to subset.
#' 
#' @return acset An acset list where every element containing feature-related
#' information has been subsetted
#'
#' @examples
#' ##load dataset
#' invisible(marinov)
#' acset = new_acset(featdata = marinov[['featdata']], refcount =
#' marinov[['refcount']], altcount = marinov[['altcount']], phenodata =
#' marinov[['phenodata']])
#'
#' ##subset features
#' sel_feat = c('HMGB1', 'ITGA4')
#' acset_filt = subset_feat(acset, sel_feat)
#' 
#' @export
subset_feat <- function(acset, feat){

    featdata = acset[['featdata']]
    
    ##subset featdata on passed feats to get passed variants
    featdata = merge(featdata, as.matrix(feat), by.x = 'feat', by.y = 1, stringsAsFactors = FALSE)
    rownames(featdata) = featdata[, 'var']

    ##subset on passed vars
    pass_var = featdata[, 'var']
    acset = subset_rows(acset, pass_var)

    return(acset)
}

#' Filter variants and features
#'
#' \code{filter_acset} removes variants and features with too little data
#'
#' The function removes variants which have less than "nmincells" cells with
#' imbalanced allelic expression. A cell is deemed to have imbalanced allelic
#' expression if its transcribed genotype is set to 0 or 2, see
#' \code{\link{call_gt}}. Subsequently it removes features with less than
#' "nminvar" variants, see \code{\link{filter_feat_nminvar}}.
#'
#' @param acset An acset list created by \code{\link{new_acset}}. It must contain
#' a "gt" element with transcribed genotype calls, see \code{\link{call_gt}}.
#' @param nmincells An integer specifying the minimum number of cells with
#' imbalanced allelic expression.
#' @param nminvar An integer specifying the minimum number of variants within a
#' feature.
#' @param feat_filter Boolean specifying if filtering should be done per feature
#' including an initial cell-filter removing cells with less than two variants
#' with imbalanced expression within the feature.
#' 
#' @return acset An acset list.
#'
#' @examples
#' ##create a small artificial genotype matrix
#' ncells = 10
#' paternal = c(0, 2, 0, 0, 2)
#' maternal = c(2, 0, 2, 2, 0)
#' gt = as.matrix(as.data.frame(rep(list(paternal, maternal), ncells / 2)))
#' vars = 1:nrow(gt)
#' colnames(gt) = 1:ncells
#' rownames(gt) = vars
#'
#' ##create a feature annotation data-frame
#' nvars = nrow(gt)
#' featdata = as.data.frame(matrix(cbind(rep('jfeat', nvars),
#' as.character(1:nvars), rep('dummy', nvars), rep('dummy', nvars)), ncol = 4,
#' dimnames = list(vars, c('feat', 'var', 'ref', 'alt'))), stringsAsFactors =
#' FALSE)
#'
#' ##create acset
#' acset = new_acset(featdata, gt = gt)
#'
#' ##Remove variants with imbalanced expression in less than nmincells cells and
#' ##features with less than nminvar variants
#' nmincells = 5
#' nminvar = 2
#' acset_filt = filter_acset(acset, nmincells, nminvar)
#' 
#' @export
filter_acset <- function(acset, nmincells = 5, nminvar = 2, feat_filter = FALSE){

    if(feat_filter){

        ##filter feats
        acset = filter_feat_gt(acset, nmincells)

    }else{ ##filter variants

        ##filter vars
        acset = filter_var_gt(acset, nmincells)

        ##filter feats on number of vars
        acset = filter_feat_nminvar(acset, nminvar)
    }

    return(acset)
}

filter_feat_counts <- function(acset, mincount = 3, ncells = 3){
    
    totcount = acset$refcount + acset$altcount
    featdata = acset$featdata

    ##loop feats
    feat2vars = tapply(featdata[, 'var'], featdata[, 'feat'], unique)
    feat_pass = unlist(lapply(feat2vars, function(jvars, totcount, mincount, ncells){
        totcount = totcount[jvars, ]
        pass_feat = filter_singlefeat_count(totcount, mincount, ncells)
        return(pass_feat)
    }, totcount, mincount, ncells))
    feat_pass = names(feat_pass)[which(feat_pass)]    
    
    acset = subset_feat(acset, feat_pass)

    return(acset)
}

filter_singlefeat_count <- function(totcount, mincount = 3, ncells = 3){
    
    ##filter cells, requring that a cell should have at least two vars with expression.
    cell2nvars = apply(totcount, 2, function(jcell, mincount){length(which(jcell >= mincount))}, mincount)
    pass_cells = names(cell2nvars)[which(cell2nvars >= 2)]

    if(length(pass_cells) >= ncells){
        totcount = totcount[, pass_cells]
        
        ##filter vars. Require at least n cells with expression
        var2ncells = apply(totcount, 1, function(jvar, mincount){length(which(jvar >= mincount))}, mincount)
        pass_var = names(var2ncells)[which(var2ncells >= ncells)]

        ##At least 2 vars
        if(length(pass_var) >= 2){
            pass_feat = TRUE
        }else{
            pass_feat = FALSE
        }
    }else{
        pass_feat = FALSE
    }
    
    return(pass_feat)
}

filter_feat_gt <- function(acset, ncells = 3){
    
    gt = acset$gt
    featdata = acset$featdata

    ##loop feats
    feat2vars = tapply(featdata[, 'var'], featdata[, 'feat'], unique)
    feat_pass = unlist(lapply(feat2vars, function(jvars, gt, ncells){
        gt = gt[jvars, ]
        pass_feat = filter_singlefeat_gt(gt, ncells)
        return(pass_feat)
    }, gt, ncells))
    feat_pass = names(feat_pass)[which(feat_pass)]    
    
    acset = subset_feat(acset, feat_pass)

    return(acset)
}

filter_singlefeat_gt <- function(gt, ncells = 3){
    
    ##filter cells, requring that a cell should have at least two vars with monoallelic calls.
    cell2nvars = apply(gt, 2, function(jcell){length(which(jcell == 0 | jcell == 2))})
    pass_cells = names(cell2nvars)[which(cell2nvars >= 2)]

    if(length(pass_cells) >= ncells){
        gt = gt[, pass_cells]
        
        ##filter vars. Require at least n cells with monoallelic calls
        var2ncells = apply(gt, 1, function(jvar){length(which(jvar == 0 | jvar == 2))})
        pass_var = names(var2ncells)[which(var2ncells >= ncells)]

        ##At least 2 vars
        if(length(pass_var) >= 2){
            pass_feat = TRUE
        }else{
            pass_feat = FALSE
        }
    }else{
        pass_feat = FALSE
    }
    
    return(pass_feat)
}

#' Filter features on minimal number of required variants
#'
#' \code{filter_feat_nminvar} removes features with less than "nmin_var"
#' variants.
#'
#' The function takes an acset and an integer, "nmin_var", as input and filter
#' out features that have less than "nmin_var" variants.
#'
#' @param acset An acset list created by the function \code{\link{new_acset}}.
#' @param nmin_var An integer specifying the minimum number of variants within
#' a feature
#' 
#' @return acset An acset list. The element [['filter']][['nmin_var']] is added
#' or updated to the acset to record the filter argument used.
#'
#' @examples
#' ##create a small artificial genotype matrix
#' ncells = 10
#' paternal = c(0, 2, 0, 0, 2)
#' maternal = c(2, 0, 2, 2, 0)
#' gt = as.matrix(as.data.frame(rep(list(paternal, maternal), ncells / 2)))
#' vars = 1:nrow(gt)
#' colnames(gt) = 1:ncells
#' rownames(gt) = vars
#'
#' ##create a feature annotation data-frame
#' nvars = nrow(gt)
#' featdata = as.data.frame(matrix(cbind(rep('jfeat', nvars),
#' as.character(1:nvars), rep('dummy', nvars), rep('dummy', nvars)), ncol = 4,
#' dimnames = list(vars, c('feat', 'var', 'ref', 'alt'))), stringsAsFactors =
#' FALSE)
#'
#' ##create acset
#' acset = new_acset(featdata, gt = gt)
#'
#' ##Remove features with less than nminvar variants
#' nminvar = 2
#' acset = filter_feat_nminvar(acset, nminvar)
#' 
#' @export
filter_feat_nminvar <- function(acset, nmin_var = 2){
###Filter on minimum number of variants present in feature
    
    ##get n variants per feature
    feat2nvar = table(acset[['featdata']][, 'feat'])

    ##get passed feats
    pass_feat = names(feat2nvar)[which(feat2nvar >= nmin_var)]

    ##subset on passed feats
    acset = subset_feat(acset, pass_feat)

    ##store filter argument
    acset[['args']][['filter']]['nmin_var'] = list(nmin_var = nmin_var)
    
    return(acset)
}

#' Filter variants
#'
#' \code{filter_var_gt} removes variants based on the number of cells with
#' imbalanced expression
#'
#' The function removes variants which have less than "nmincells" cells with
#' imbalanced allelic expression. A cell is deemed to have imbalanced allelic
#' expression if its transcribed genotype is set to 0 or 2, see
#' \code{\link{call_gt}}.
#'
#' @param acset An acset list created by \code{\link{new_acset}}. It must
#' contain a "gt" element with transcribed genotype calls, see
#' \code{\link{call_gt}}.
#' @param nmincells An integer specifying the minimum number of cells with
#' imbalanced allelic expression.
#' 
#' @return acset An acset list subsetted on variants that pass the filter.
#'
#' @examples
#' ##create a small artificial genotype matrix
#' ncells = 10
#' paternal = c(0, 2, 0, 0, 2)
#' maternal = c(2, 0, 2, 2, 0)
#' gt = as.matrix(as.data.frame(rep(list(paternal, maternal), ncells / 2)))
#' vars = 1:nrow(gt)
#' colnames(gt) = 1:ncells
#' rownames(gt) = vars
#'
#' ##create a feature annotation data-frame
#' nvars = nrow(gt)
#' featdata = as.data.frame(matrix(cbind(rep('jfeat', nvars),
#' as.character(1:nvars), rep('dummy', nvars), rep('dummy', nvars)), ncol = 4,
#' dimnames = list(vars, c('feat', 'var', 'ref', 'alt'))), stringsAsFactors =
#' FALSE)
#'
#' ##create acset
#' acset = new_acset(featdata, gt = gt)
#'
#' ##Remove variants with imbalanced expression in less than 3 cells
#' nmincells = 3
#' acset_filt = filter_var_gt(acset, nmincells)
#' 
#' @export
filter_var_gt <- function(acset, nmincells = 3){
    
    gt = acset$gt
    featdata = acset$featdata

    ##loop feats
    feat2vars = tapply(featdata[, 'var'], featdata[, 'feat'], unique)
    var_pass = unlist(lapply(feat2vars, function(jvars, gt, nmincells){
        gt = gt[jvars, , drop = FALSE]
        pass_var = filter_var_singlefeat_gt(gt, nmincells)
        return(pass_var)
    }, gt, nmincells))
    
    acset = subset_rows(acset, var_pass)

    return(acset)
}

filter_var_singlefeat_gt <- function(gt, ncells = 3){
    
    ##filter cells, requring that a cell should have at least two vars with monoallelic calls.
    cell2nvars = apply(gt, 2, function(jcell){length(which(jcell == 0 | jcell == 2))})
    pass_cells = names(cell2nvars)[which(cell2nvars >= 2)]

    if(length(pass_cells) >= ncells){
        gt = gt[, pass_cells]
        
        ##filter vars. Require at least n cells with monoallelic calls
        var2ncells = apply(gt, 1, function(jvar){length(which(jvar == 0 | jvar == 2))})
        pass_var = names(var2ncells)[which(var2ncells >= ncells)]

    }else{
        pass_var = NULL
    }
    
    return(pass_var)
}

filter_var_mincount <- function(acset, mincount = 3, ncells = 1){
##Filter vars on expression in at least ncells    
    totcount = acset$refcount + acset$altcount
    var2ncells = apply(totcount, 1, function(jvar, mincount){length(which(jvar >= mincount))}, mincount)
    pass_vars = names(var2ncells)[which(var2ncells >= ncells)]
    acset = subset_rows(acset, pass_vars)

    return(acset)
}

#' Filter homozygous variants
#'
#' \code{filter_homovars} removes variants for which a large proportion of cells
#' with imbalanced expression is imbalanced towards the same allele
#'
#' The function removes variants which tend to express the same allele in a
#' large proportion of cells, that is, the allelic expression is stable across
#' cells rather than random with respect to which allele is the most highly
#' expressed. The purpose of this filter is to reduce the number of variants
#' that have falsely been called as heterozygous variants but are actually
#' homozygous. Use with caution as a heterozygous variants can indeed have
#' imbalanced expression towards the same allele in the majority of cells.
#' A cell is deemed to have imbalanced allelic expression if its allele-specific
#' expression, ase, is < mono_ase or >(1 - mono_ase), where ase = alternative
#' allele count / (alternative allele count + reference allele count). If the
#' number of cells expressing one allele in such an imbalanced manner is
#' significantly greater than the number of cells expressing the other allele
#' the variant is removed (binomial test).
#'
#' @param acset An acset list which must contain "refcount" and "altcount"
#' elements with allele counts, see \code{\link{new_acset}}.
#' @param alpha A numeric specifying the significance level of the binomial
#' test, where the test compares the number of cells expressing each allele
#' in an imbalanced "monoallelic" manner.
#' @param mono_ase A numeric between 0 and 1 specifying the allele specific
#' expression level at which to deem an allele monoallelically expressed.
#' 
#' @return acset An acset list subsetted on variants that pass the filter.
#'
#' @examples
#' ##load dataset
#' invisible(marinov)
#' acset = new_acset(featdata = marinov[['featdata']], refcount =
#' marinov[['refcount']], altcount = marinov[['altcount']], phenodata =
#' marinov[['phenodata']])
#'
#' ##Remove variants having monoallelic expression of the same allele in a
#' ##large proportion of cells
#' alpha = 0.1
#' mono_ase = 0.1
#' acset_filt = filter_homovars(acset, alpha, mono_ase)
#' 
#' @export
filter_homovars <- function(acset, alpha = 0.1, mono_ase = 0.1){

    ##Number of cells with ase towards each of the alleles
    ase = acset[['altcount']] / (acset[['refcount']] + acset[['altcount']])
    n.refcells = apply(ase, 1, function(j.ase){length(which(j.ase < mono_ase))})
    n.altcells = apply(ase, 1, function(j.ase){length(which(j.ase > (1 - mono_ase)))})
    n.allele2cells = cbind(n.refcells, n.altcells, n.refcells + n.altcells)
    colnames(n.allele2cells) = c('ref', 'alt', 'sum')

    bi.rs = rownames(n.allele2cells)[which(n.allele2cells[, 'sum'] == 0)]
    n.allele2cells = n.allele2cells[setdiff(rownames(n.allele2cells), bi.rs), ]
    
    ##get binom pval
    pvals = apply(n.allele2cells, 1, function(j.c){pvals = stats::binom.test(j.c['alt'], j.c['sum'], p = 0.5, alternative = 'two.sided')$p.value})
    pass.rs = names(pvals)[which(pvals >= alpha)]
    pass.rs = c(pass.rs, bi.rs)
    
    ##filter
    acset = subset_rows(acset, pass.rs);
    
    return(acset)
}

#' Filter variants with zero allele counts in all cells
#'
#' \code{filter_zerorow} removes variants which do no have any read count
#'
#' @param acset An acset list which must contain "refcount" and "altcount"
#' elements with allele counts, see \code{\link{new_acset}}.
#' 
#' @return acset An acset list subsetted on variants that pass the filter.
#'
#' @examples
#' ##create small artifical allele count matrixes where first row is zero
#' ncells = 10
#' ref = c(0, 10, 1, 5, 2, 9)
#' alt = c(0, 2, 5, 5, 6, 3)        
#' refcount = as.matrix(as.data.frame(rep(list(ref, alt), ncells / 2)))
#' altcount = as.matrix(as.data.frame(rep(list(alt, ref), ncells / 2)))
#' 
#' vars = 1:nrow(refcount)
#' samples = 1:ncells
#' colnames(refcount) = samples
#' rownames(refcount) = vars
#' colnames(altcount) = samples
#' rownames(altcount) = vars
#'
#' ##create a feature annotation data-frame
#' nvars = length(vars)
#' featdata = as.data.frame(matrix(cbind(rep('jfeat', nvars),
#' as.character(1:nvars), rep('dummy', nvars), rep('dummy', nvars)), ncol = 4,
#' dimnames = list(vars, c('feat', 'var', 'ref', 'alt'))), stringsAsFactors =
#' FALSE)
#'
#' ##create acset
#' acset = new_acset(featdata, altcount = altcount, refcount = refcount)
#'
#' ##Remove variants having no allele counts
#' acset_filt = filter_zerorow(acset)
#' 
#' @export
filter_zerorow <- function(acset){

    totcount = acset[['altcount']] + acset[['refcount']]
    rowcount = apply(totcount, 1, sum)
    pass_rs = names(rowcount)[which(rowcount != 0)]
    acset = subset_rows(acset, pass_rs)

    return(acset)
}

#' Call transcribed genotypes
#'
#' \code{call_gt} calls transcribed genotypes using allele RNA-seq read counts
#'
#' This is a simplistic transcribed genotype caller which is used to discretize
#' which allele is the most expressed. For each variant and cell one of four
#' possible discrete values is set depending on the expression of the two
#' alleles. 0: reference allele most highly expressed, 1: bi-allelic expression
#' with similar degree of expression from the two alleles, 2: alternative allele
#' most highly expressed, or, NA if there are not enough reads to make a call.
#'
#' @param acset An acset list created by the function \code{\link{new_acset}}.
#' The acset must contain a refcount and altcount matrix with allele counts.
#' @param min_acount An integer specifying the minimum number of reads required
#' to be present for at least one of the two alleles as to make a call. If not
#' fulfilled the call is set to NA.
#' @param fc An integer specifying the fold-change cutoff between the expression
#' of the alternative and reference allele (fc.observed = alternative allele
#' count / reference allele count). The transcribed genotype call is set to 0 if
#' fc.observed <= fc, 2 if fc.observed >= fc and 1 otherwise (if there are
#' enough reads present, see the "min_acount" parameter).
#' 
#' @return acset An acset list where the "gt" element is set or updated. "gt" is
#' a matrix with integer values representing transcribed genotype calls. 0:
#' reference allele most highly expressed, 1: bi-allelic expression with similar
#' degree of expression from the two alleles, 2: alternative allele most highly
#' expressed. NA's are used to represent entries where no call could be made.
#' The rownames are set to the variant column, "var", of featdata. The colnames
#' are set to the colnames of "refcount". The elements
#' [['args']][['filter']][c('min_acount', 'fc')] are added or updated to the
#' acset as to record the filter arguments used.
#'
#' @examples
#' ##load dataset
#' invisible(marinov)
#' acset = new_acset(featdata = marinov[['featdata']], refcount =
#' marinov[['refcount']], altcount = marinov[['altcount']], phenodata =
#' marinov[['phenodata']])
#'
#' ##Call transcribed genotypes
#' min_acount = 3
#' fc = 3
#' acset = call_gt(acset, min_acount, fc)
#' 
#' @export
call_gt <- function(acset, min_acount = 3, fc = 3){
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
    samples = colnames(acset[['refcount']])
    gt = matrix(data = NA, nrow = length(feats), ncol = length(samples), dimnames = list(feats, samples))
    gt[r.mono.ind] = 0
    gt[bi.ind] = 1
    gt[a.mono.ind] = 2
    
    acset[['gt']] = gt

    ##set complementary genotype matrix
    acset = set_compl_gt(acset)

    ##store filter argument
    acset[['args']][['filter']][c('min_acount', 'fc')] = list(min_acount = min_acount, fc = fc)

    return(acset)
}

set_haplotype <- function(acset){
##NB: The phasing is performend independently per feature, therefore, between two features we do not if there are switches in the haplotype or not.

    featdata = acset[['featdata']]
    varflip = acset[['varflip']]
    
    vars = rownames(featdata)
    hapA = featdata[, 'ref']
    hapB = featdata[, 'alt']
    names(hapA) = vars
    names(hapB) = vars
    hapA[varflip] = featdata[varflip, 'alt']
    hapB[varflip] = featdata[varflip, 'ref']
    
    phasedfeat = featdata[, c('feat', 'var', 'ref', 'alt')]
    phasedfeat = cbind(phasedfeat, hapA, hapB, stringsAsFactors = FALSE)

    acset[['phasedfeat']] = phasedfeat
    
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
