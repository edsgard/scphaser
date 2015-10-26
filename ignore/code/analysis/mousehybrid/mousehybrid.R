

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

    ##genotype
    min_acount = 10
    fc = 4
    acset = call_gt(acset, min_acount, fc)

    ##filter vars
    nmincells = 10
    acset = filter_var_gt(acset, nmincells)
    lapply(acset, dim) #2370

    ##filter feats
    nminvar = 2
    acset = filter_feat_nminvar(acset, nminvar)
    lapply(acset, dim) #2343
    length(unique(acset$featdata$feat)) #299

    
    ##*###
    ##Performance
    ##*###

    ##Create synthetic data where truth is known
    nfracflip = 0.5
    input = 'gt' ##'ac'
    weigh = FALSE ##TRUE
    method = 'exhaust' ##'cluster'
    synt_res = synt_flip(acset$featdata, acset$gt, nfracflip, input, weigh, method)

    ##performance
    varsflip = synt_res[['varsflip']]
    perf = get_performance(varsflip, acset$featdata)

    ##MCC: 
    ##FPR: 1.5%
    ##TPR: 97.6%
    ##nfeats_id: 282 / 299 genes
    
}

subset_chr <- function(acset, chr){
    featdata = acset$featdata
    pass_vars = featdata[which(featdata$chr == chr), 'var']
    acset_sub = subset_rows(acset, pass_vars)

    return(acset_sub)
}

synt_flip <- function(featdata, gt, nfracflip, input, weigh, method){
##start with gt where all vars are initially phased (no vars should be flipped)
##then randomly flip nfracflip variants
    
    ##Randomly flip nfracflip vars
    vars = featdata$var
    nvars = length(vars)
    nflip = round(nfracflip * nvars)
    vars_flip = sample(vars, nflip)

    gt_compl = compl_gt(gt)
    gt[vars_flip, ] = gt_compl[vars_flip, ]

    acset_synt = new_acset(featdata = featdata, gt = gt)

    ##phase
    acset_synt = phase(acset_synt, input = input, weigh = weigh, method = method)

    
    ##*###
    ##Set observed flips to the complement if closer to exp_compl
    ##*###

    ##binarize vars known to be flipped (expected)
    vars_flip_bin = as.integer(vars %in% vars_flip)
    varsflip_exp = cbind(vars, featdata[vars, 'feat'], vars_flip_bin)
    colnames(varsflip_exp) = c('var', 'feat', 'exp')
    varsflip = as.data.frame(varsflip_exp, stringsAsFactors = FALSE)
    rownames(varsflip) = varsflip[, 'var']
    varsflip[, 'exp'] = as.integer(varsflip[, 'exp'])
    
    ##add complement of expected flip
    exp_compl = varsflip[, 'exp']
    exp_compl[varsflip[, 'exp'] == 0] = 1
    exp_compl[varsflip[, 'exp'] == 1] = 0
    varsflip = cbind(varsflip, exp_compl)
    
    ##subset
    varsflip = varsflip[acset_synt$featdata$var, ]

    ##add column of observed flips
    obs_vars = acset_synt[['varflip']]
    obs = as.integer(varsflip[, 'var'] %in% obs_vars)
    varsflip = cbind(varsflip, obs)

    ##set the adjusted observed flip-vector to the complementary of obs, if obs is closer to exp_compl
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

    
    return(list(varsflip = varsflip, acset = acset_synt))
}

get_performance <- function(varsflip, featdata){
    
    n_true = sum(varsflip[, 'exp'])
    
    varsflip_neg = varsflip[which(varsflip[, 'obs_adj'] == 0), ]
    n_neg = nrow(varsflip_neg)
    tn = length(which(varsflip_neg[, 'exp'] == 0))
    fn = n_neg - tn
    
    varsflip_pos = varsflip[which(varsflip[, 'obs_adj'] == 1), ]
    n_pos = sum(varsflip[, 'obs_adj'])
    tp = sum(varsflip_pos[, 'exp'])
    fp = n_pos - tp
    
    fpr = fp / n_pos
    tpr = tp / n_true
    mcc = ((tp * tn) - (fp * fn)) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    
    ##check at feat-level, may be better.
    feat2vars = tapply(featdata[, 'var'], featdata[, 'feat'], unique)
    varsflip[, 'obs_adj'] = as.integer(varsflip[, 'obs_adj'])
    
    feat2correct = unlist(lapply(feat2vars, function(vars, varsflip){
        jflip = varsflip[vars, ]
        identical(jflip[, 'exp'], jflip[, 'obs_adj'])
    }, varsflip))
    n_id = length(which(feat2correct)) #342
    nfeats = length(feat2vars) #691

    var_perf = c(mcc, fpr, tpr, n_true, n_pos, tp, fp, tn, fn)
    names(var_perf) = c('mcc', 'fpr', 'tpr', 'true', 'pos', 'tp', 'fp', 'tn', 'fn')
    feat_perf = c(nfeats, n_id)
    names(feat_perf) = c('n_feats', 'n_identical')
    
    perf = list(var = var_perf, feat = feat_perf)
    return(perf)
}
