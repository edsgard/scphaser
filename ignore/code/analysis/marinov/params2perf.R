
###Assess complexity under varying parameters
###NB: The paths are relative to the root of the scphaser git repo
###
###SYNOPSIS
###library('devtools')
###devtools::load_all()
###source("ignore/code/analysis/marinov/params2perf.R")


##Params
sys = 'dna'


##*###
##Dirs
##*###

if(sys == 'dna'){
    cloud.dir = '/mnt/kauffman/edsgard/cloud/btsync/work/rspd'
}
if(sys == 'rs13'){
    cloud.dir = '/Volumes/Data/cloud/btsync/work/rspd'
}

##OUT
out_rds_dir = './ignore/res/marinov/data'
out_pdf_dir = './ignore/res/marinov/pdf'


##*###
##Files
##*###

##IN
ref.counts.rds = file.path(cloud.dir, 'projects/scphaser/nogit/data/marinov', paste('ref', '.counts.rds', sep = ''))
alt.counts.rds = file.path(cloud.dir, 'projects/scphaser/nogit/data/marinov', paste('alt', '.counts.rds', sep = ''))
snp.annot.rds = file.path(cloud.dir, 'projects/scphaser/nogit/data/marinov', 'snp.annot.rds')

##OUT
phased_rds = file.path(out_rds_dir, 'phased.acset.rds')
times.pdf = file.path(out_pdf_dir, 'phasetime.bar.pdf')

##Libs
source('./ignore/code/analysis/performance.R')
library('BiocParallel')
library('dplyr')
library('tidyr')

main <- function(){


    ##*###
    ##Read and prep data
    ##*###
    altcount = readRDS(alt.counts.rds)
    refcount = readRDS(ref.counts.rds)
    featdata = readRDS(snp.annot.rds)

    var = featdata[, 'rsid']
    feat = featdata[, 'gene']
    featdata = cbind(featdata, var, feat, stringsAsFactors = FALSE)
    
    ##create acset
    acset = new_acset(featdata, refcount, altcount)
    lapply(acset, dim) #43,313    28
    length(unique(acset[['featdata']][, 'gene'])) #9,456
    feat2nvars = table(acset[['featdata']][, 'feat'])
    table(feat2nvars) #16: 51

    
    ##*###
    ##Params
    ##*###
    
    ##Call gt
    min_acount = c(3, 5, 10)
    fc = c(3, 4, 5) #75/25
    
    ##Filter variants on n.cells monoallelic and feats with < 2 j.vars
    nmincells = c(3, 5, 10)
    nminvar = 2 #default
    
    ##phasing method params
    input = c('gt', 'ac')
    method = c('exhaust', 'cluster')
    weigh = c(TRUE, FALSE)    

    ##specify paramset as all possible combinations of the params
    paramset = expand.grid(min_acount, fc, nmincells, nminvar, input, weigh, method, stringsAsFactors = FALSE)
    colnames(paramset) = c('min_acount', 'fc', 'nmincells', 'nminvar', 'input', 'weigh', 'method')

    ##parallelization
    ncores = 80
    bp_param = BiocParallel::MulticoreParam(workers = ncores)

    
    ##*###
    ##Phase
    ##*###
    ##loop param set
    nparamset = nrow(paramset)
    phased_list = lapply(1:nparamset, filter_phase_par, paramset = paramset, acset = acset, bp_param = bp_param)

    ##Dump
    saveRDS(phased_list, file = phased_rds)
    ##status: sub (12.36 -> 12.36; 80 cores). #data: 3,739x28; paramset: 8 rows.

    
    
}
