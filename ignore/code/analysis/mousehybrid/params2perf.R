
###Assess performance under varying parameters
###NB: The paths are relative to the root of the scphaser git repo
###
###SYNOPSIS
###library('devtools')
###devtools::load_all()
###source("ignore/code/analysis/mousehybrid/params2perf.R")


##Dirs
out_data_dir = './ignore/res/mousehybrid/perf/data'
out_pdf_dir = './ignore/res/mousehybrid/perf/pdf'

##Files
perf_rds = file.path(out_data_dir, 'perf.full.rds')

##Libs
source('./ignore/code/analysis/performance.R')
library('BiocParallel')
library('dplyr')
library('tidyr')

main <- function(){


    ##*###
    ##Read and prep data
    ##*###
    load('../nogit/data/mousehybrid/mousehybrid.rda')

    ##create acset
    acset = new_acset(mousehybrid$featdata, mousehybrid$refcount, mousehybrid$altcount, mousehybrid$phenodata)
    lapply(acset, dim)
    
    ##subset
    ngenes = 100
    acset = subset_ngenes(acset, ngenes)

    
    ##*###
    ##Specify params
    ##*###

    ##*##
    ##Small test
    ##*##

    verbosity = -1
    
    ##phasing method params
    input = 'gt' #c('gt', 'ac')
    method = c('exhaust', 'cluster')
    weigh = FALSE #c(TRUE, FALSE)    
        
    ##permutation iterations
    npermiter = 2 #10
    perm_iter = 1:npermiter

    ##min read count for genotyping
    min_acount = c(3, 5, 10)
    
    ##allelic fold-change for mono-allelic call
    fc = 4 # c(3, 4, 5)
    
    ##filter vars on min number of cells with mono-allelic call
    nmincells = 10 #c(3, 5, 10)

    ##specify paramset as all possible combinations of the params
    paramset = expand.grid(perm_iter, min_acount, fc, nmincells, input, weigh, method, stringsAsFactors = FALSE)
    colnames(paramset) = c('perm_iter', 'min_acount', 'fc', 'nmincells', 'input', 'weigh', 'method')

    ##parallelization
    ncores = 2
    bp_param = BiocParallel::MulticoreParam(workers = ncores)


    ##*###
    ##Full run
    ##*###

    verbosity = 0
    
    ##phasing method params
    input = c('gt', 'ac')
    method = c('exhaust', 'cluster')
    weigh = c(TRUE, FALSE)    
        
    ##permutation iterations
    npermiter = 10
    perm_iter = 1:npermiter

    ##min read count for genotyping
    min_acount = c(3, 5, 10)
    
    ##allelic fold-change for mono-allelic call
    fc = c(3, 4, 5)
    
    ##filter vars on min number of cells with mono-allelic call
    nmincells = c(3, 5, 10)

    ##specify paramset as all possible combinations of the params
    paramset = expand.grid(perm_iter, min_acount, fc, nmincells, input, weigh, method, stringsAsFactors = FALSE)
    colnames(paramset) = c('perm_iter', 'min_acount', 'fc', 'nmincells', 'input', 'weigh', 'method')

    ##parallelization
    ncores = 80
    bp_param = BiocParallel::MulticoreParam(workers = ncores)

    
    
    ##*###
    ##Calculate performance
    ##*###
    ##loop phase_params and get a dataframe of performance with respect to the other params for each phasing parameter set

    ##loop param set
    nparamset = nrow(paramset)
    perf_list = BiocParallel::bplapply(1:nparamset, filter_get_perf_par, BPPARAM = bp_param, paramset = paramset, acset = acset, verbosity = verbosity)
    perf_df = do.call(rbind, perf_list)

    ##bind params and perf
    params2perf_df = cbind(paramset, perf_df)    

    ##*###
    ##Number of retained vars and feats wrt filter parameters
    ##*###
    filter2nvars_nfeats = get_postfilter_stats(acset, paramset)
    params2perf_df = cbind(params2perf_df, filter2nvars_nfeats)
    params2perf_df = tidyr::unite(params2perf_df, in.meth.w, input, method, weigh, remove = FALSE, sep = '.')

    ##Dump
    saveRDS(params2perf_df, file = perf_rds)
    ##status: fin (2016.02.10: 15.38->16.12 (34min); 80 cores). #paramset: 2160 rows.
    
        
    
    ##*###
    ##Plot performance
    ##*###

    ##read data
    params2perf_df = readRDS(perf_rds)

    ##summarize by method, across filters
    a = dplyr::group_by(params2perf_df, input, method, weigh)
    dplyr::summarize(a, mcc = mean(mcc))

    ##add numeric 
    weigh.num = as.numeric(params2perf_df[, 'weigh'])
    params2perf_df = cbind(params2perf_df, weigh.num)
    which(unlist(lapply(params2perf_df, is.factor)))
        
    obs.long = params2perf_df
    x.str = 'min_acount'
    group.col = 'in.meth.w'
    colour.col = 'method'
    lt.col = 'input'
    size.col = 'weigh.num'
    j.pdf = file.path(out_pdf_dir, 'mcc.pdf')
    pdf(j.pdf)
    plot_perf(obs.long, x.str = x, group.col = group.col, colour.col = colour.col, lt.col = lt.col, size.col = size.col)
    dev.off()


    ##*###
    ##TPR vs FPR
    ##*###

    ##select one method (the winner)
    sel.meth = 'ac.exhaust.TRUE'
    obs.long = dplyr::filter(params2perf_df, in.meth.w == sel.meth)

    ##num -> factors
    obs.long[, 'nmincells'] = factor(obs.long[, 'nmincells'])
    obs.long[, 'min_acount'] = factor(obs.long[, 'min_acount'])
        
    ##barplot
    pdf.h = 6
    pdf.w = 10
    j.pdf = file.path(out_pdf_dir, sel.meth, 'tpr.fpr.bar.pdf')
    dir.create(dirname(j.pdf))
    pdf(j.pdf, width = pdf.w, height = pdf.h)
    plot_tprfpr(obs.long)
    dev.off()

    
    ##Obsolete: mean fpr and tpr
    obs.long = obs.long %>% group_by(min_acount, fc, nmincells) %>% summarise(tpr.mean = mean(tpr), fpr.mean = mean(fpr), tpr.sd = sd(tpr), fpr.sd = sd(fpr))

    x.str = 'fpr.mean'
    y.str = 'tpr.mean'
    colour.col = 'nmincells'
    group.col = 'nmincells'
    shape.col = 'min_acount'
    
    gg = ggplot(obs.long, aes_string(x = x.str, y = y.str, colour = colour.col, group = group.col, shape = shape.col))
    gg = gg + geom_point()
    gg = gg + geom_line()
    gg = gg + facet_grid(~fc)
    print(gg)
      

    ##*###
    ##Filters vs "sensitivity"
    ##*###
    frac.vars = obs.long[, 'nvars_postfilt'] / obs.long[, 'nvars_prefilt']
    obs.long = cbind(obs.long, frac.vars)
    fracvars2filters = unique(obs.long[, c('fc', 'nmincells', 'min_acount', 'frac.vars')])
    fracvars2filters = tidyr::unite(fracvars2filters, fc.minacount, fc, min_acount, remove = FALSE, sep = '.')
    fracvars2filters[, 'min_acount'] = factor(fracvars2filters[, 'min_acount'])
    fracvars2filters[, 'fc'] = factor(fracvars2filters[, 'fc'])
    
    x.str = 'nmincells'
    y.str = 'frac.vars'
    colour.col = 'fc'
    lt.col = 'min_acount'
    group.col = 'fc.minacount'
    y.lab = 'Fraction variants'
    txt.size = 15
    j.pdf = file.path(out_pdf_dir, 'filtparam2sens.pdf')
    pdf(j.pdf)
    plot_filtsens(fracvars2filters, x.str, y.str, colour.col, lt.col, group.col, y.lab = y.lab, txt.size = txt.size)
    dev.off()
    

    ##*###
    ##Number of errors (as in marinov/post.phase.R)
    ##*###
    library('ggplot2')
    n.err = params2perf_df[, 'fp'] + params2perf_df[, 'fn']
    frac.err = n.err / params2perf_df[, 'nvars_postfilt']
    params2perf_df = cbind(params2perf_df, n.err, frac.err)

    y.resp = c('Number of vars', 'Fraction erronously phased vars', 'Number of genes')
    names(y.resp) = c('nvars_postfilt', 'frac.err', 'nfeats_postfilt')
    
    for(j.y.resp.name in names(y.resp)){
        param2perf = params2perf_df
        
        ##add numeric weigh col
        weigh.num = as.numeric(param2perf[, 'weigh'])
        param2perf = cbind(param2perf, weigh.num)
        which(unlist(lapply(param2perf, is.factor)))

        x.str = 'min_acount'
        group.col = 'in.meth.w'
        colour.col = 'method'
        lt.col = 'input'
        size.col = 'weigh.num'
        obs.long = param2perf
        x.lab = 'Minimum allele count of at least one allele'
        y.lab = paste(y.resp[j.y.resp.name], sep = ' ')
        
        ##data and mappings
        gg = ggplot(obs.long, aes_string(x = x.str, y = j.y.resp.name, colour = colour.col, group = group.col, linetype = lt.col, size = size.col))
        gg = gg + scale_size(range = c(0.5, 1), breaks = c(0, 1), labels = c('FALSE', 'TRUE'))

        ##median value as line
        gg = gg + stat_summary(fun.y = median, geom = "line")

        ##95% CI
        gg = gg + stat_summary(data = obs.long, fun.data = mean_cl_boot, geom = "linerange") ## geom = 'errorbar'

        ##facet
        gg = gg + facet_grid(fc ~ nmincells, scales = 'free_y')
        
        ##Background
        gg = gg + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(panel.background = element_blank())

        ##x ticks
        gg = gg + theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"), axis.ticks = element_line(colour = 'black'))

        ##gg = gg + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
        
        ##axis lables
        gg = gg + xlab(x.lab) + ylab(y.lab)

        j.pdf = file.path(out_pdf_dir, paste(j.y.resp.name, 'pdf', sep = '.'))
        pdf(file = j.pdf)
        print(gg)
        dev.off()
    }    
    
    
    ##*##
    ##TBD: matrix of plots for each phasing param.set
    ##*###
    ##split into list
    perf_list = split(params2perf_df, params2perf_df[, 'in.meth.w'])

    paramsets = names(perf_list)
    for(jset in paramsets){

        jperf = perf_list[[jset]]
        
        j.pdf = file.path(pdf.dir, paste(jset, '.mcc.pdf', sep = ''))
        pdf(file = j.pdf)
        plot_perf(jperf)
        dev.off()
    }

}
