
###Assess performance under varying parameters
###NB: The paths are relative to the root of the scphaser git repo
###
###SYNOPSIS
###library('devtools')
###devtools::load_all()
###source("ignore/code/analysis/mousehybrid/params2perf.R")


##Dirs
out_data_dir = './ignore/res/perf/data'
out_pdf_dir = './ignore/res/perf/pdf'

##Files
perf_rds = file.path(out_data_dir, 'perf.rds')

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
    
    ##phasing method params
    input = 'gt' #c('gt', 'ac')
    method = c('exhaust', 'cluster')
    weigh = FALSE #c(TRUE, FALSE)    
        
    ##permutation iterations
    npermiter = 2
    perm_iter = 1:npermiter

    ##min read count for genotyping
    min_acount = 10 #c(3, 10)
    
    ##allelic fold-change for mono-allelic call
    fc = 4 # c(3, 4)
    
    ##filter vars on min number of cells with mono-allelic call
    nmincells = 10 #c(5, 10)

    ##specify paramset as all possible combinations of the params
    paramset = expand.grid(perm_iter, min_acount, fc, nmincells, input, weigh, method, stringsAsFactors = FALSE)
    colnames(paramset) = c('perm_iter', 'min_acount', 'fc', 'nmincells', 'input', 'weigh', 'method')

    ##parallelization
    ncores = 2
    bp_param = BiocParallel::MulticoreParam(workers = ncores)

    
    ##*###
    ##Calculate performance
    ##*###
    ##loop phase_params and get a dataframe of performance with respect to the other params for each phasing parameter set

    ##loop param set
    nparamset = nrow(paramset)
    perf_list = BiocParallel::bplapply(1:nparamset, filter_get_perf_par, BPPARAM = bp_param, paramset = paramset, acset = acset)
    perf_df = do.call(rbind, perf_list)

    ##bind params and perf
    params2perf_df = cbind(paramset, perf_df)    
    
    ##Dump
    saveRDS(params2perf_df, file = perf_rds)

    
    ##*###
    ##Number of retained vars and feats wrt filter parameters
    ##*###
    filter2nvars_nfeats = get_postfilter_stats(acset, paramset)
    params2perf_df = cbind(params2perf_df, filter2nvars_nfeats)
    
    
    ##*###
    ##Plot performance
    ##*###
    ##TODO: matrix of plots for each phasing parameter set

    params2perf_df = readRDS(perf_rds)

    ##summarize by method, across filters
    a = dplyr::group_by(params2perf_df, input, method, weigh)
    dplyr::summarize(a, mcc = mean(mcc))

    ##split into list
    params2perf_df = tidyr::unite(params2perf_df, in.meth.w, input, method, weigh, remove = FALSE, sep = '.')
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
