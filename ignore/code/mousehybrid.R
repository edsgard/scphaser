

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

    
    ##*###
    ##Subset
    ##*###
    featdata = acset$featdata
    pass_vars = featdata[which(featdata$chr == 'chr1'), 'var']
    acset_sub = subset_rows(acset, pass_vars)
        
    ##phase
    acset_sub = phase(acset_sub, input = 'gt', weigh = FALSE, method = 'exhaust')

    ##nothing should be flipped:
    dim(acset_sub$featdata) #11521
    length(acset_sub$varflip) #139 vars

    
    ##*###
    ##Full set
    ##*###
    acset = phase(acset, input = 'gt', weigh = FALSE, method = 'exhaust')
    
    ##FPR at variant level
    featdata = acset$featdata
    nvars = length(featdata$var) ##167443
    nvar_flipped = length(acset$varflip) ##2034
    nvar_flipped / nvars ##1.2%

    ##FPR at feat level
    featdata = acset$featdata
    nfeats_flipped = length(unique(featdata[acset$varflip, 'feat']))
    nfeats = length(unique(featdata$feat))
    nfeats_flipped / nfeats ##5.7%
    
    ##gene score before phasing
    score_prephasing = get_var_gt_mat(acset_sub)

    
    ##*###
    ##FPR and TPR
    ##*###    
    ##flip some variants (true cases) and check if refinds them

    ##Synth data
    featdata = acset$featdata
    vars = featdata$var
    vars_flip = sample(vars, length(vars) / 2)

    gt = acset$gt
    gt_compl = acset$gt_compl
    gt[vars_flip, ] = gt_compl[vars_flip, ]

    acset_synt = new_acset(featdata = featdata, gt = gt)

    ##subset
    jchr = 'chr1'
    acset_sub = subset_chr(acset_synt, jchr)

    ##phase
    acset_sub = phase(acset_sub, input = 'gt', weigh = FALSE, method = 'exhaust')

    ##compare flipped vars with those that should have been flipped
    jchr_vars = featdata[which(featdata[, 'chr'] == jchr), 'var']


    ##*###
    ##Set observed flips to the flip-set closest to exp or exp_compl
    ##*###
    vars_flip_bin = as.integer(vars %in% vars_flip)
    varsflip_exp = cbind(vars, featdata[vars, 'feat'], vars_flip_bin)
    colnames(varsflip_exp) = c('var', 'feat', 'exp')
    varsflip = as.data.frame(varsflip_exp, stringsAsFactors = FALSE)
    rownames(varsflip) = varsflip[, 'var']
    varsflip[, 'exp'] = as.integer(varsflip[, 'exp'])
    
    ##add compl
    exp_compl = varsflip[, 'exp']
    exp_compl[varsflip[, 'exp'] == 0] = 1
    exp_compl[varsflip[, 'exp'] == 1] = 0
    varsflip = cbind(varsflip, exp_compl)
    
    ##subset
    varsflip = varsflip[acset_sub$featdata$var, ]

    ##add column of observed flips
    obs_vars = acset_sub[['varflip']]
    obs = as.integer(varsflip[, 'var'] %in% obs_vars)
    varsflip = cbind(varsflip, obs)

    ##hamming distance of obs to both exp and exp_compl, within each feat
    feat2vars = tapply(varsflip[, 'var'], varsflip[, 'feat'], unique)
    varsflip_list = lapply(feat2vars, function(vars, varsflip){
        jflip = varsflip[vars, ]
        exp_error = length(which(jflip[, 'obs'] != jflip[, 'exp']))
        exp_compl_error = length(which(jflip[, 'obs'] != jflip[, 'exp_compl']))
        obs_adj = jflip[, 'obs']
        if(exp_compl_error < exp_error){
            obs = jflip[, 'obs']
            obs_adj[which(obs == 0)] = 1
            obs_adj[which(obs == 1)] = 0
        }
        jflip = cbind(jflip, obs_adj)

        return(jflip)
    }, varsflip)
    varsflip = do.call(rbind, varsflip_list)
    rownames(varsflip) = varsflip[, 'var']

    
    n_true = sum(varsflip[, 'exp']) #5719
    n_pos = sum(varsflip[, 'obs_adj']) #5238
    varsflip_pos = varsflip[which(varsflip[, 'obs_adj'] == 1), ]
    n_tp = sum(varsflip_pos[, 'exp'])
    n_fp = n_pos - n_tp
    fpr = n_fp / n_pos #19.4%
    tpr = n_tp / n_true #73.8%

    ##check at feat-level, may be better.
    featsub = acset_sub$featdata
    feat2vars = tapply(featsub[, 'var'], featsub[, 'feat'], unique)
    varsflip[, 'obs_adj'] = as.integer(varsflip[, 'obs_adj'])
    
    feat2correct = unlist(lapply(feat2vars, function(vars, varsflip){
        jflip = varsflip[vars, ]
        identical(jflip[, 'exp'], jflip[, 'obs_adj'])
    }, varsflip))
    n_id = length(which(feat2correct)) #342
    nfeats = length(unique(featsub[, 'feat'])) #691
    
}

subset_chr <- function(acset, chr){
    featdata = acset$featdata
    pass_vars = featdata[which(featdata$chr == chr), 'var']
    acset_sub = subset_rows(acset, pass_vars)

    return(acset_sub)
}
