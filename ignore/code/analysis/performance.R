

plot_perf <- function(obs.long, x.str = 'min_acount', group.col = 'in.meth.w', colour.col = 'method', lt.col = 'input', size.col = 'weigh.num', y.lab = 'MCC'){

    library('ggplot2')
    
    ##data and mappings
    gg = ggplot(obs.long, aes_string(x = x.str, y = 'mcc', colour = colour.col, group = group.col, linetype = lt.col, size = size.col))

    gg = gg + scale_size(range = c(0.5, 1), breaks = c(0, 1), labels = c('FALSE', 'TRUE'))
    
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
    ##gg = gg + ylab(y.lab)

    print(gg)
    ##colors
    ##gg = gg + scale_color_manual(values = lin2color)

}

plot_tprfpr <- function(obs.long, errbar.w = 0.25){

    fill.col = 'nmincells'
    dodge = position_dodge(width = 0.9)

    obs.long[, 'fpr'] = -obs.long[, 'fpr']

    gg = ggplot(obs.long, aes_string(x = 'min_acount', y = 'tpr', fill = fill.col))
    gg = gg + geom_bar(stat = 'summary', fun.y = 'mean', position = dodge)
    gg = gg + geom_errorbar(stat = 'summary', fun.data = mean_se, position = dodge, width = errbar.w)
    gg = gg + geom_bar(aes_string(y = 'fpr'), stat = 'summary', fun.y = 'mean', position = dodge)
    gg = gg + geom_errorbar(aes_string(y = 'fpr'), stat = 'summary', fun.data = mean_se, position = dodge, width = errbar.w)
    gg = gg + facet_grid(~fc)
    ##gg = gg + coord_cartesian(ylim = c(-0.05, 1))
    gg = gg + scale_y_continuous(breaks = c(-0.05, seq(0, 1, 0.25)))
    gg = gg + geom_hline(yintercept = 0)
    gg = gg + geom_hline(yintercept = 0.95, linetype = 'dashed', colour = 'grey')
    gg = gg + geom_hline(yintercept = -0.05, linetype = 'dashed', colour = 'grey')
    
    ##Background
    gg = gg + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(panel.background = element_blank())
    gg = gg + theme(axis.text = element_text(colour="black"), axis.ticks = element_line(colour = 'black'))

    print(gg)
}

plot_filtsens <- function(obs.long, x.str = 'nmincells', y.str = 'frac.vars', colour.col = 'fc', lt.col = 'min_acount', group.col = 'fc.minacount', y.lab = 'Fraction variants', txt.size = 10){

    library('ggplot2')
    
    ##data and mappings
    gg = ggplot(obs.long, aes_string(x = x.str, y = y.str, colour = colour.col, group = group.col, linetype = lt.col))
    gg = gg + geom_line()
    
    ##facet
    ##gg = gg + facet_grid(fc ~ nmincells, scales = 'free_y')
    
    ##Background
    gg = gg + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(panel.background = element_blank())

    ##x ticks
    gg = gg + theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"), axis.ticks = element_line(colour = 'black'))

    gg = gg + theme(text = element_text(size = txt.size))

    ##axis limits
    gg = gg + coord_cartesian(ylim = c(0, 1))
    
    ##axis lables
    gg = gg + ylab(y.lab)

    print(gg)
    ##colors
    ##gg = gg + scale_color_manual(values = lin2color)

}

get_perf_paramset <- function(acset, perm_iter, min_acount, fc, nmincells, input, weigh, method, feat_filter = FALSE, bp_param = BiocParallel::SerialParam()){
    
    ##specify paramset as all possible combinations of the params
    paramset = expand.grid(perm_iter, min_acount, fc, nmincells, stringsAsFactors = FALSE)
    colnames(paramset) = c('perm_iter', 'min_acount', 'fc', 'nmincells')
    nparamset = nrow(paramset)

    ##loop param set
    perf_list_var = list()
    perf_list_feat = list()
    length(perf_list_var) = nparamset
    length(perf_list_feat) = nparamset
    for(jparam in 1:nparamset){
        jparamset = paramset[jparam, ]
        jperf = filter_get_perf(acset, jparamset[, 'min_acount'], jparamset[, 'fc'], jparamset[, 'nmincells'], input, weigh, method, feat_filter = feat_filter, bp_param = bp_param)
        perf_list_var[[jparam]] = jperf[['var']]
        perf_list_feat[[jparam]] = jperf[['feat']]
    }
    perf_df_var = do.call(rbind, perf_list_var)
    perf_df_feat = do.call(rbind, perf_list_feat)
    
    ##bind params and perf
    params2perf_df_var = cbind(paramset, perf_df_var)
    params2perf_df_feat = cbind(paramset, perf_df_feat)
    param2perf_list = list(var = params2perf_df_var, feat = params2perf_df_feat)
    
    return(params2perf_df_list)
}

filter_get_perf_par <- function(jparam, paramset, acset, input, weigh, method, feat_filter = FALSE, bp_param = BiocParallel::SerialParam(), verbosity = 0){

    verbose = R.utils::Verbose(threshold = -1)
    R.utils::cat(verbose, jparam)
    jparamset = paramset[jparam, ]
    jperf = filter_get_perf(acset, jparamset[, 'min_acount'], jparamset[, 'fc'], jparamset[, 'nmincells'], jparamset[, 'input'], jparamset[, 'weigh'], jparamset[, 'method'], feat_filter = feat_filter, bp_param = bp_param, verbosity = verbosity)
        
    return(jperf)    
}

filter_phase_par <- function(jparam, paramset, acset, feat_filter = FALSE, bp_param = BiocParallel::SerialParam(), verbosity = 0){

    verbose = R.utils::Verbose(threshold = -1)
    R.utils::cat(verbose, jparam)
    jparamset = paramset[jparam, ]


    ##Call gt
    min_acount = jparamset[, 'min_acount']
    fc = jparamset[, 'fc']
    acset = call_gt(acset, min_acount, fc)

    ##Filter
    nmincells = jparamset[, 'nmincells']
    nminvar = jparamset[, 'nminvar']
    acset = filter_acset(acset, nmincells, nminvar)
    
    ##phase
    input = jparamset[, 'input']
    weigh = jparamset[, 'weigh']
    method = jparamset[, 'method']

    start.time = proc.time()
    acset_phased = phase(acset, input = input, weigh = weigh, method = method, bp_param = bp_param, verbosity = verbosity)
    end.time = proc.time()
    pass.time = end.time - start.time

    ##bind time and call params
    acset_phased[['time']] = pass.time
    acset_phased[['params']] = jparamset
    
    return(acset_phased)    
}

phase_par <- function(jparam, paramset, acset, feat_filter = FALSE, bp_param = BiocParallel::SerialParam(), verbosity = 0){

    verbose = R.utils::Verbose(threshold = -1)
    R.utils::cat(verbose, jparam)
    jparamset = paramset[jparam, ]

    input = jparamset[, 'input']
    weigh = jparamset[, 'weigh']
    method = jparamset[, 'method']

    ##phase
    start.time = proc.time()
    acset_phased = phase(acset, input = input, weigh = weigh, method = method, bp_param = bp_param, verbosity = verbosity)
    end.time = proc.time()
    pass.time = end.time - start.time

    ##bind time and call params
    acset_phased[['time']] = pass.time
    acset_phased[['params']] = jparamset
    
    return(acset_phased)    
}

complexity_par <- function(jparam, paramset, acset, bp_param = BiocParallel::SerialParam(), verbosity = 0){

    ##print progress
    verbose = R.utils::Verbose(threshold = -1)
    R.utils::cat(verbose, jparam)
    jparamset = paramset[jparam, ]

    ##*###
    ##filter
    ##*###
        
    ##subsample n.cells
    j.cells = jparamset[, 'n.cells']
    cells = acset[['phenodata']][, 'sample']
    pass_cells = sample(cells, j.cells)
    acset_filt = subset_cols(acset, pass_cells)

    ##Filter variants on n.cells monoallelic and feats with < 2 j.vars
    nmincells = jparamset[, 'nmincells']
    j.vars = jparamset[, 'n.vars']    
    acset_filt = filter_acset(acset_filt, nmincells, j.vars)

    ##Filter feats on != j.vars (TBD: could also subsample)
    featdata = acset_filt[['featdata']]
    feat2nvars = table(featdata[, 'feat'])
    pass_feats = names(feat2nvars)[which(feat2nvars == j.vars)]
    
    ##subsample n.sel feats
    j.feats = jparamset[, 'n.feats']
    pass_feats = sample(pass_feats, j.feats, replace = FALSE)
    acset_filt = subset_feat(acset_filt, pass_feats)

    ##get dim post-filtering
    filt.dim = dim(acset_filt[['refcount']])
    names(filt.dim) = c('n.row', 'n.col')

    
    ##*###
    ##phase
    ##*###
    input = jparamset[, 'input']
    weigh = jparamset[, 'weigh']
    method = jparamset[, 'method']
    start.time = proc.time()
    acset_dummy = phase(acset_filt, input = input, weigh = weigh, method = method, bp_param = bp_param, verbosity = verbosity)
    end.time = proc.time()
    pass.time = end.time - start.time
    
    ##add time to params
    param.time = c(jparamset, filt.dim, pass.time)
    
    return(param.time)
}

get_postfilter_stats <- function(acset, paramset){
    
    nparamset = nrow(paramset)
    featdata = acset$featdata
    nvars = nrow(featdata)
    nfeats = length(unique(featdata[, 'feat']))
    nvars_prefilt = rep(nvars, nparamset)
    nfeats_prefilt = rep(nfeats, nparamset)
    nvars_postfilt = rep(NA, nparamset)
    nfeats_postfilt = rep(NA, nparamset)
    for(jparam_it in 1:nparamset){
        print(jparam_it)
        jparamset = paramset[jparam_it, ]
        
        ##Call gt
        acset_filt = call_gt(acset, jparamset[, 'min_acount'], jparamset[, 'fc'])

        ##Filter
        acset_filt = filter_acset(acset_filt, jparamset[, 'nmincells'])

        featdata = acset_filt[['featdata']]
        nvars_postfilt[jparam_it]= nrow(featdata)
        nfeats_postfilt[jparam_it] = length(unique(featdata[, 'feat']))        
    }

    filter2nvars = cbind(nfeats_prefilt, nvars_prefilt, nfeats_postfilt, nvars_postfilt)
    
    return(filter2nvars)
}

filter_get_perf <- function(acset, min_acount, fc, nmincells, input, weigh, method, nminvar = 2, nfracflip = 0.5, feat_filter = FALSE, bp_param = BiocParallel::SerialParam(), verbosity = 0){
    
    ##Call gt
    acset = call_gt(acset, min_acount, fc)

    ##Filter
    acset = filter_acset(acset, nmincells, nminvar, feat_filter)
            
    ##Create synthetic data where truth is known
    synt_res = synt_flip(acset, nfracflip, input, weigh, method, bp_param, verbosity = verbosity)

    ##performance
    varsflip = synt_res[['varsflip']]
    perf = get_performance(varsflip, acset$featdata)

    return(perf)
}

rand_flip <- function(acset, nfracflip = 0.5){
    
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

    return(acset_synt)
}

synt_flip <- function(acset, nfracflip, input, weigh, method, bp_param = BiocParallel::SerialParam(), verbosity = 0){
##start with gt where all vars are initially phased (no vars should be flipped)
##then randomly flip nfracflip variants

    ##randomly flip nfracflip vars
    acset_synt = rand_flip(acset, nfracflip)
        
    ##phase
    acset_synt = phase(acset_synt, input = input, weigh = weigh, method = method, bp_param = bp_param, verbosity = verbosity)

    
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
