
###Assess complexity under varying parameters
###NB: The paths are relative to the root of the scphaser git repo
###
###SYNOPSIS
###library('devtools')
###devtools::load_all()
###source("ignore/code/analysis/mousehybrid/complexity.R")


##Dirs
out_data_dir = './ignore/res/complexity/data'
out_pdf_dir = './ignore/res/complexity/pdf'

##Files
timing_rds = file.path(out_data_dir, 'timing.rds')

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
    lapply(acset, dim) #167443    336
    length(unique(acset[['featdata']][, 'feat'])) #11590
    feat2nvars = table(acset[['featdata']][, 'feat'])
    table(feat2nvars) #16: 333

    
    ##*###
    ##Parameters
    ##*###

    ##*###
    ##Test run
    ##*###
    verbosity = 0
    
    ##sampling iterations
    perm_iter = 1
    
    ##phasing method params
    input = c('gt')
    method = c('exhaust', 'cluster')
    weigh = c(FALSE)
    
    ##"minimal" filtering
    min_acount = 3
    fc = 3
    nmincells = 3
    nminvar = 2
    
    ##filter on n.vars
    n.vars = c(2)

    ##subsample n.sel feats
    n.feats = c(25)
    
    ##subsample n.cells
    n.cells = c(50)    

    ##specify paramset as all possible combinations of the params
    paramset = expand.grid(perm_iter, min_acount, fc, nmincells, input, weigh, method, n.vars, n.feats, n.cells, stringsAsFactors = FALSE)
    colnames(paramset) = c('perm_iter', 'min_acount', 'fc', 'nmincells', 'input', 'weigh', 'method', 'n.vars', 'n.feats', 'n.cells')
    dim(paramset) #2
    
    ##parallelization
    ncores = 2
    bp_param = BiocParallel::MulticoreParam(workers = ncores)

    
    ##*###
    ##Full run
    ##*###
    verbosity = 0
    
    ##sampling iterations
    perm_iter = 1
    
    ##phasing method params
    input = c('gt', 'ac')
    method = c('exhaust', 'cluster')
    weigh = c(TRUE, FALSE)
    
    ##"minimal" filtering
    min_acount = 3
    fc = 3
    nmincells = 3
    nminvar = 2
    
    ##filter on n.vars
    n.vars = c(2, 4, 8, 16)

    ##subsample n.sel feats
    n.feats = c(25, 50, 100)
    
    ##subsample n.cells
    n.cells = c(50, 100, 200, 300)    

    ##specify paramset as all possible combinations of the params
    paramset = expand.grid(perm_iter, min_acount, fc, nmincells, input, weigh, method, n.vars, n.feats, n.cells, stringsAsFactors = FALSE)
    colnames(paramset) = c('perm_iter', 'min_acount', 'fc', 'nmincells', 'input', 'weigh', 'method', 'n.vars', 'n.feats', 'n.cells')
    dim(paramset) #384
    
    ##parallelization
    ncores = 80
    bp_param = BiocParallel::MulticoreParam(workers = ncores)


    ##*###
    ##Filter
    ##*###
        
    ##Call gt
    acset = call_gt(acset, min_acount, fc)

    ##Filter
    acset = filter_acset(acset, nmincells, nminvar)
    lapply(acset, dim) #167443 -> 152110
    length(unique(acset[['featdata']][, 'feat'])) #11590
    
    feat2nvars = table(acset[['featdata']][, 'feat'])
    table(feat2nvars) #16: 294    
        
    
    ##*###
    ##phase (parallelize on the paramset)
    ##*###
    ##time each phasing

    ##loop param set
    nparamset = nrow(paramset)
    res_list = BiocParallel::bplapply(1:nparamset, complexity_par, BPPARAM = bp_param, paramset = paramset, acset = acset, verbosity = verbosity)
    res_df = do.call(rbind, res_list)
    res_df = as.data.frame(res_df)

    ##convert datatypes
    char.cols = c('input', 'weigh', 'method')
    num.cols = c('user.self', 'sys.self', 'elapsed', 'user.child', 'sys.child')
    int.cols = setdiff(colnames(res_df), c(char.cols, num.cols))
    res_df[, char.cols] = lapply(res_df[, char.cols], as.character)
    res_df[, num.cols] = lapply(res_df[, num.cols], as.numeric)
    res_df[, int.cols] = lapply(res_df[, int.cols], as.integer)
    
    ##add combo column
    res_df = tidyr::unite(res_df, in.meth.w, input, method, weigh, remove = FALSE, sep = '.')

    ##Dump
    saveRDS(res_df, file = timing_rds)
    ##status: sub (11.05 -> ; 384 rows, 80 cores)
    ##Error in `row.names<-.data.frame`(`*tmp*`, value = value): duplicate 'row.names' are not allowed
    
    
    ##*###
    ##Post-phase analysis (plotting)
    ##*###    
    res_df = readRDS(timing_rds)

    x.str = 'n.cells'
    y.str = 'elapsed'
    group.col = 'in.meth.w'
    colour.col = 'method'
    lt.col = 'input'
    size.col = 'weigh.num'
    
    gg = ggplot(res_df, aes_string(x = x.str, y = y.str, group = group.col, colour = colour.col, lt = lt.col, size = size.col))
    gg = gg + geom_line()
    gg = gg + scale_size(range = c(0.5, 1), breaks = c(0, 1), labels = c('FALSE', 'TRUE'))    
    gg = gg + facet_grid(n.vars ~ n.feats, scales = 'free_y')


    
}

obsolete <- function(){

    ##test filtering
    feat2nvars = as.data.frame(feat2nvars, stringsAsFactors = FALSE)
    colnames(feat2nvars) = c('feat', 'n.vars')    
    feat2nvars = merge(feat2nvars, as.matrix(sel.nvars), by.x = 'n.vars', by.y = 1)
    lapply(feat2nvars, class)
    feat.subset = feat2nvars %>% group_by(n.vars) %>% sample_n(n.sel)

}
