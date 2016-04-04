
###SYNOPSIS
###library('devtools')
###devtools::load_all()
###source("ignore/code/analysis/humantcell/phase.R")

##Dirs
out_rds_dir = './ignore/res/humantcell/rds'
out_pdf_dir = './ignore/res/humantcell/pdf'

##Files
phased_rds = file.path(out_rds_dir, 'phased.acset.rds')
times.pdf = file.path(out_pdf_dir, 'phasetime.bar.pdf')

##Libs
source('./ignore/code/analysis/performance.R')
library('reshape2')
library('tidyr')
library('ggplot2')

main <- function(){


    ##*###
    ##Read and prep data
    ##*###
    load('data/humantcell.rda')

    ##create acset
    acset = new_acset(humantcell$featdata, humantcell$refcount, humantcell$altcount, humantcell$phenodata)
    lapply(acset, dim) #872, 505

    ##Call gt
    min_acount = c(3, 5, 10)
    fc = c(3, 4, 5) #75/25
    acset = call_gt(acset, min_acount, fc)
    
    ##Filter variants on n.cells monoallelic and feats with < 2 j.vars
    nmincells = c(3, 5, 10)
    nminvar = 2    
    acset = filter_acset(acset, nmincells, nminvar)
    lapply(acset, dim) #802, 505
    length(unique(acset[['featdata']][, 'feat'])) #292
    
    ##phasing method params
    input = c('gt', 'ac')
    method = c('exhaust', 'cluster')
    weigh = c(TRUE, FALSE)    

    ##specify paramset as all possible combinations of the params
    paramset = expand.grid(min_acount, fc, nmincells, input, weigh, method, stringsAsFactors = FALSE)
    colnames(paramset) = c('min_acount', 'fc', 'nmincells', 'input', 'weigh', 'method')
    
    ##parallelization
    ncores = 80
    bp_param = BiocParallel::MulticoreParam(workers = ncores)

    
    ##*###
    ##Phase
    ##*###
    ##loop param set
    nparamset = nrow(paramset)
    phased_list = lapply(1:nparamset, phase_par, paramset = paramset, acset = acset, bp_param = bp_param)
    
    ##Dump
    saveRDS(phased_list, file = phased_rds)
    ##status: sub (16.16->16.17; 80 cores). #data: 802x505; paramset: 8 rows.

    
    ##*###
    ##Post-phase analysis
    ##*###
    phased.list = readRDS(phased_rds)

    lapply(phased.list, function(j.acset){lapply(j.acset, dim)})

    ##get params and time
    params2time = lapply(phased.list, function(j.acset){c(j.acset$time, j.acset$params)})
    params2time.df = as.data.frame(do.call(rbind, params2time))
    num.cols = c('user.self', 'sys.self', 'elapsed', 'user.child', 'sys.child')
    int.cols = c('min_acount', 'fc', 'nmincells')
    char.cols = c('input', 'weigh', 'method')
    params2time.df[, num.cols] = lapply(params2time.df[, num.cols], as.numeric)
    params2time.df[, int.cols] = lapply(params2time.df[, int.cols], as.integer)
    params2time.df[, char.cols] = lapply(params2time.df[, char.cols], as.character)
    params2time.df = tidyr::unite(params2time.df, in.meth.w, input, method, weigh, remove = FALSE, sep = '.')

    ##order by time
    params2time.df = params2time.df[order(params2time.df[, 'elapsed']), ]
    params2time.df[, 'in.meth.w'] = factor(params2time.df[, 'in.meth.w'], levels = params2time.df[, 'in.meth.w'])
    
    ##plot
    x.str = 'in.meth.w'
    y.str = 'elapsed'
    gg = ggplot(params2time.df, aes_string(x = x.str, y = y.str))
    gg = gg + geom_bar(stat = 'identity')
    gg = gg + ylab('Time') + xlab('')
    
    ##Background
    gg = gg + theme(axis.text.x = element_text(angle = 90), axis.text = element_text(size = 20), axis.title = element_text(size = 24))
    gg = gg + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(panel.background = element_blank())
    gg = gg + theme(axis.text = element_text(colour="black"), axis.ticks = element_line(colour = 'black'))

    pdf(times.pdf)
    plot(gg)
    dev.off()
    
    
}
