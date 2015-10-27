

plot_perf <- function(obs.long, your.x, your.colour, y.lab = 'MCC'){

    ##data and mappings
    gg = ggplot(obs.long, aes_string(x = 'min_acount', y = 'mcc', colour = your.colour, group = group.col))

    ##facet
    gg = gg + facet_grid(fc ~ nmincells, scales = 'free_y')
    
    ##mean value as line
    gg = gg + stat_summary(fun.y = median, geom = "line")

    ##95% CI
    gg = gg + stat_summary(data = obs.long, fun.data = mean_cl_boot, geom = "linerange") ## geom = 'errorbar'
    
    ##Background
    gg = gg + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(panel.background = element_blank())

    ##x ticks
    gg = gg + theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"), axis.ticks = element_line(colour = 'black'))

    ##gg = gg + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    
    ##axis lables
    gg = gg + ylab(y.lab) ##+ xlab(x.lab)

    ##colors
    ##gg = gg + scale_color_manual(values = lin2color)

}

get_perf_paramset <- function(acset, perm_iter, min_acount, fc, nmincells, input, weigh, method){

    
    ##specify paramset as all possible combinations of the params
    paramset = expand.grid(perm_iter, min_acount, fc, nmincells, stringsAsFactors = FALSE)
    colnames(paramset) = c('perm_iter', 'min_acount', 'fc', 'nmincells')
    nparamset = nrow(paramset)

    ##loop param set
    perf_list = list()
    length(perf_list) = nparamset
    for(jparam in 1:nparamset){
        jparamset = paramset[jparam, ]
        jperf = filter_get_perf(acset, jparamset[, 'min_acount'], jparamset[, 'fc'], jparamset[, 'nmincells'], input, weigh, method)
        perf_list[[jparam]] = jperf[['var']]
    }
    perf_df = do.call(rbind, perf_list)

    ##bind params and perf
    params2perf_df = cbind(paramset, perf_df)
        
    return(params2perf_df)
}

filter_get_perf <- function(acset, min_acount, fc, nmincells, input, weigh, method, nminvar = 2, nfracflip = 0.5, feat_filter = TRUE){


    ##*###
    ##genotype
    ##*###
    acset = call_gt(acset, min_acount, fc)


    ##*###
    ##Filter
    ##*###

    if(feat_filter){

        ##filter feats on 
        acset = filter_feat_gt(acset, nmincells)

    }else{ ##filter variants

        ##filter feats on 
        acset = filter_var_gt(acset, nmincells)

        ##filter feats
        acset = filter_feat_nminvar(acset, nminvar)
    }


    
    ##*###
    ##Performance
    ##*###
    
    ##Create synthetic data where truth is known
    synt_res = synt_flip(acset, nfracflip, input, weigh, method)

    ##performance
    varsflip = synt_res[['varsflip']]
    perf = get_performance(varsflip, acset$featdata)

    return(perf)
}

synt_flip <- function(acset, nfracflip, input, weigh, method){
##start with gt where all vars are initially phased (no vars should be flipped)
##then randomly flip nfracflip variants

    featdata = acset$featdata
    gt = acset$gt
    
    ##Randomly flip nfracflip vars
    vars = featdata$var
    nvars = length(vars)
    nflip = round(nfracflip * nvars)
    vars_flip = sample(vars, nflip)

    gt_compl = compl_gt(gt)
    gt[vars_flip, ] = gt_compl[vars_flip, ]
    
    if(!is.null(acset$refcount)){
        refcount = acset$refcount
        altcount = acset$altcount
        newref = refcount
        newalt = altcount
        newref[vars_flip, ] = altcount[vars_flip, ]
        newalt[vars_flip, ] = refcount[vars_flip, ]
        
        acset_synt = new_acset(featdata = featdata, gt = gt, refcount = newref, altcount = newalt)
    }else{
        acset_synt = new_acset(featdata = featdata, gt = gt)
    }
    
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

subset_ngenes <- function(acset, ngenes){
    
    featdata = acset$featdata
    feat2var = tapply(featdata$var, featdata$feat, unique)
    sel_vars = unlist(feat2var[1:ngenes])
    acset = subset_rows(acset, sel_vars)

    return(acset)
}
subset_chr <- function(acset, chr){
    featdata = acset$featdata
    pass_vars = featdata[which(featdata$chr == chr), 'var']
    acset_sub = subset_rows(acset, pass_vars)

    return(acset_sub)
}
