


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

score_var <- function(acset){
###for each variant calculate the "distance" to a phased reference variant

    gt = acset[['gt_phased']]
    ##TODO: check if gt_phased exists

    ##exclude columns with only NAs and/or 1s    
    n_na = apply(gt, 2, function(jcell){length(which(is.na(jcell) | jcell == 1))})
    ncells = ncol(gt)
    keep_ind = which(n_na != ncells) 
    gt = gt[, keep_ind]
    
    ##create reference variant from the most frequent gt within each cell
    ref_var = apply(gt, 2, get_mostfreqgt)

    ##get agreement btw vars and ref_var
    nshared = apply(gt, 1, get_ngtshared, ref_var = ref_var)
    fracshared = nshared / ncol(gt)

    fracshared = t(fracshared)
    
    return(fracshared)
}

get_ngtshared <- function(jvar, ref_var){
    
    na_ind = which(is.na(jvar))
    bi_ind = which(jvar == 1)
    gt_ind = setdiff(1:length(jvar), na_ind)

    shared = rep(0, 3)
    names(shared) = c('gt', 'na', 'bi')
    
    shared['gt'] = length(which(jvar[gt_ind] == ref_var[gt_ind]))
    shared['na'] = length(na_ind)
    shared['bi'] = length(bi_ind)

    return(shared)
}

get_mostfreqgt <- function(jcell){
###Return the most frequent item, after removal of NAs.
###Also, assume that 1's should be excluded.
    
    jcell = jcell[!(is.na(jcell) | jcell == 1)]
    
    gt2nvars = sort(table(jcell), decreasing = TRUE)
    mostfreq = names(gt2nvars)[1];

    return(mostfreq)
}

get_phase_pval <- function(acset, nperm = 100, fixedrowmargin = FALSE){

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
            acset_rnd = call_gt(acset_rnd, min_acount, fc)
        }else{
            acset_rnd = racset(acset, type = 'gt', fixedrowmargin = fixedrowmargin)
        }
                
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

#' Randomize an acset
#'
#' \code{racset} randomizes the entries in the allele count matrixes or the
#' genotype matrix
#'
#' The function scrambles all entries within each allele count matrix or within
#' the genotype matrix.
#'
#' @param acset An acset list created by \code{\link{new_acset}}.
#' @param type A character string with two allowed values. 'ac' specifies that
#' the allele count matrixes should be randomized and 'gt' specifies that the
#' genotype matrix should be randomized.
#' @param fixedrowmargin Boolean specifying if the row-sums should be kept
#' unchanged.
#' 
#' @return acset A randomized acset list.
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
#' ##Randomize the genotype matrix
#' type = 'gt'
#' acset_rand = racset(acset, type)
#' 
#' @export
racset <- function(acset, type = 'ac', fixedrowmargin = FALSE){

    if(type == 'ac'){
        acset_rnd = rac(acset, fixedrowmargin)
    }
    if(type == 'gt'){
        acset_rnd = rgt(acset, fixedrowmargin)
    }    
    
    return(acset_rnd)
}

rac <- function(acset, fixedrowmargin = FALSE){

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

rgt <- function(acset, fixedrowmargin = FALSE){
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

#' Phasing concordance plot
#'
#' \code{plot_conc} plots the phasing concordance, see \code{\link{set_gt_conc}}
#'
#' To assess the success of the phasing one can calculate the degree of
#' variability remaining if all cells with haplotype 2 are set to haplotype 1.
#' As a rough measure of this this function calculates the variability as the
#' number of cells that differ from the inferred haplotype for every gene with
#' two variants. The differing number of cells per gene we denote as the
#' inconcordance and the number of cells with identical haplotype to the
#' inferred haplotype as concordance, see \code{\link{set_gt_conc}} for further
#' details.
#'
#' @param acset An acset list created by \code{\link{new_acset}}. The acset must
#' contain the elements 'gt_conc' and 'gt_phased_conc' which contain the
#' concordance and inconcordance for every gene before (gt_conc) and after
#' phasing (gt_phased_conc), see \code{\link{set_gt_conc}}.
#' @param feats A character vector of feature names to include in the plot.
#' @param cex A numerical value giving the amount by which plotting text and
#' symbols should be magnified relative to the default, see \code{\link{par}}.
#'
#' @examples
#' ##load dataset
#' invisible(marinov)
#' acset = new_acset(featdata = marinov[['featdata']], refcount =
#' marinov[['refcount']], altcount = marinov[['altcount']], phenodata =
#' marinov[['phenodata']])
#'
#' ##Call gt
#' acset = call_gt(acset, min_acount = 3, fc = 3)
#'
#' ##Filter variants and genes
#' acset = filter_acset(acset, nmincells = 5, nminvar = 2)
#'
#' ##Phase
#' acset = phase(acset, input = 'gt', weigh = FALSE, method = 'exhaust')
#'
#' ##Calculate the genotype concordance before and after phasing
#' acset = set_gt_conc(acset)
#' 
#' ##Plot the genotype concordance before and after phasing
#' acset = plot_conc(acset)
#' 
#' @export
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
    x_txt = c(20, 20)
    
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

#' Cell-distribution variability upon phasing
#'
#' \code{set_gt_conc} calculates an approximate measure of the spread of the
#' cells in the allele-specific expression variant-space per gene
#'
#' To assess the success of the phasing one can calculate the degree of
#' variability remaining if all cells with haplotype 2 are set to haplotype 1.
#' As a rough measure of this this function calculates the variability as the
#' number of cells that differ from the inferred haplotype for every gene with
#' two variants. The differing number of cells per gene we denote as the
#' inconcordance and the number of cells with identical haplotype to the
#' inferred haplotype as concordance. This can be understood by viewing each
#' cell as a point in a transcribed genotype variant space, where each dimension
#' is a variant (within a gene), with values in the transcribed genotype domain.
#' The expression of the alleles within a gene from a single cell will then tend
#' to cluster towards the haplotype state, and the remaining variability of the
#' cells in that space after phasing can be used to measure how well the cells
#' conform to the haplotype.
#'
#' @param acset An acset list created by \code{\link{new_acset}}. The acset must
#' contain an element 'gt', see \code{\link{call_gt}} and an element
#' 'gt_phased', see function \code{\link{phase}}.
#' 
#' @return acset An acset list where two elements have been added or updated,
#' 'gt_conc' and 'gt_phased_conc'. The elements contain the concordance and
#' inconcordance for every gene before (gt_conc) and after phasing
#' (gt_phased_conc).
#'
#' @examples
#' ##load dataset
#' invisible(marinov)
#' acset = new_acset(featdata = marinov[['featdata']], refcount =
#' marinov[['refcount']], altcount = marinov[['altcount']], phenodata =
#' marinov[['phenodata']])
#'
#' ##Call gt
#' acset = call_gt(acset, min_acount = 3, fc = 3)
#'
#' ##Filter variants and genes
#' acset = filter_acset(acset, nmincells = 5, nminvar = 2)
#'
#' ##Phase
#' acset = phase(acset, input = 'gt', weigh = FALSE, method = 'exhaust')
#'
#' ##Get genotype concordance before and after phasing
#' acset = set_gt_conc(acset)
#' 
#' @export
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
