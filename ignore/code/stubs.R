
get_initial_gt <- function(acset){
###We want to start with the matrix which induce an increase in the compactness score
###TBD: Don't think I'm actually doing something sensible here! Because, if the number of ref monoallelic calls are greater why can't we just as well start with that?!

    ##not finished...
    
    ##Calculate the fraction of alternative monoallelic calls, for each SNP
    n_samples = ncol(jfeat_gt)
    var2fracalt = rowSums(jfeat_gt == 2) / n_samples

    ##Calculate the fraction of SNPs with var2fracalt >0.5
    n_vars = nrow(jfeat_gt)
    frac_vars = length(which(var2fracalt > 0.5)) / n_vars
    
    ##If frac_vars > 0.5 start out with SNPs with var2fracalt < 0.5, otherwise the other set of SNPs.
    if(frac_vars > 0.5){
        jfeat_gt_phased = jfeat_gt
    }else{
        jfeat_gt_phased = jfeat_gt_comple
    }

    return(jfeat_gt_phased)
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
