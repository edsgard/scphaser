

##NB: The paths are relative to the root of the scphaser git repo

source('./ignore/code/analysis/performance.R')

main <- function(){

    ##Read data
    load('../nogit/data/mousehybrid/mousehybrid.rda')

    ##create acset
    acset = new_acset(mousehybrid$featdata, mousehybrid$refcount, mousehybrid$altcount, mousehybrid$phenodata)
    lapply(acset, dim)

    ##genotype
    min_acount = 3
    fc = 3
    acset = call_gt(acset, min_acount, fc)

    
    ##*######################
    ##Full set
    ##*######################
    acset = phase(acset, input = 'gt', weigh = FALSE, method = 'exhaust')
    
    ##Flipped variants
    featdata = acset$featdata
    nvars = length(featdata$var) ##167443
    nvar_flipped = length(acset$varflip) ##2034
    nvar_flipped / nvars #1.2%

    ##Feats with flipped variants
    featdata = acset$featdata
    nfeats_flipped = length(unique(featdata[acset$varflip, 'feat']))
    nfeats = length(unique(featdata$feat))
    nfeats_flipped / nfeats ##5.7%
    
    ##gene score before phasing
    score_prephasing = get_var_gt_mat(acset_sub)

    
    ##*######################
    ##Subset
    ##*######################
    featdata = acset$featdata
    pass_vars = featdata[which(featdata$chr == 'chr1'), 'var']
    acset_sub = subset_rows(acset, pass_vars)

    ##genotype
    min_acount = 3
    fc = 3
    acset_sub = call_gt(acset_sub, min_acount, fc)

    ##phase
    acset_sub = phase(acset_sub, input = 'gt', weigh = FALSE, method = 'exhaust')

    ##nothing should be flipped:
    nvars = length(acset_sub$featdata$var) #11521
    nvar_flipped = length(acset_sub$varflip) #139 vars
    nvar_flipped / nvars #1.2%    

    ##*###
    ##Performance
    ##*###
    ##flip some variants (true cases) and check if we can phase them

    ##Create synthetic data where truth is known
    nfracflip = 0.5
    synt_res = synt_flip(acset_sub$featdata, acset_sub$gt, nfracflip)
    
    ##Get performance
    varsflip = synt_res[['varsflip']]
    perf = get_performance(varsflip, acset_sub$featdata)
    
    ##FPR: 21%
    ##TPR: 75%
    ##nfeats_id: 342 / 691 genes

    ##difference to consensus among cells where not NA or bi.
    varscores = score_var(acset)
    maxfrac = 1 - apply(varscores[, c('na', 'bi')], 1, sum)
    fracofconsensus = varscores[, 'gt'] / maxfrac
    summary(fracofconsensus)

    

    ##*#########################
    ##Harder filters
    ##*#########################
    
    ##create acset
    acset = new_acset(mousehybrid$featdata, mousehybrid$refcount, mousehybrid$altcount, mousehybrid$phenodata)
    lapply(acset, dim)

    ##subset
    jchr = 'chr1'
    acset = subset_chr(acset, jchr)
    lapply(acset, dim) #11521
    length(unique(acset$featdata$feat)) #691

    min_acount = 10
    fc = 4
    nmincells = 10
    input = 'gt' ##'ac'
    weigh = FALSE ##TRUE
    method = 'exhaust' ##'cluster'

    perf = filter_get_perf(acset, min_acount, fc, nmincells, input, weigh, method)
    print(perf)
    ##MCC: 
    ##FPR: 1.5%
    ##TPR: 97.6%
    ##nfeats_id: 282 / 299 genes
            
}

