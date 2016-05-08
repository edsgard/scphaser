
###Assess complexity under varying parameters
###NB: The paths are relative to the root of the scphaser git repo
###
###SYNOPSIS
###library('devtools')
###devtools::load_all()
###source("ignore/code/analysis/borel/phase.R")


##Params
sys = 'dna'
sys = 'rs13'

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
out_rds_dir = '../nogit/data/borel'
out_pdf_dir = './ignore/res/borel/pdf'


##*###
##Files
##*###

##IN
ref.counts.rds = file.path(cloud.dir, 'projects/scphaser/nogit/data/borel', paste('ref', '.counts.rds', sep = ''))
alt.counts.rds = file.path(cloud.dir, 'projects/scphaser/nogit/data/borel', paste('alt', '.counts.rds', sep = ''))
snp.annot.rds = file.path(cloud.dir, 'projects/scphaser/nogit/data/borel', 'snp.annot.rds')

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
    
    ##create acset
    acset = new_acset(featdata, refcount, altcount)
    lapply(acset, dim) #    163
    length(unique(acset[['featdata']][, 'feat'])) #13,698
    feat2nvars = table(acset[['featdata']][, 'feat'])
    table(feat2nvars) #16: 159

    ##Filter vars with 0 counts
    acset = filter_zerorow(acset)
    lapply(acset, dim) #230,674    163
    
    ##Call gt
    min_acount = 3
    fc = 3 #75/25
    acset = call_gt(acset, min_acount, fc)

    ##Filter variants on number of cells where ASE towards the same allele
    alpha = 0.1
    mono.ase = 0.1
    if(!(mono.ase == 0)){
        acset = filter_homovars(acset, alpha = alpha, mono.ase = mono.ase)
    }
    lapply(acset, dim) #220,793    163

    ##Filter variants on n.cells monoallelic and feats with < 2 j.vars
    nmincells = 5
    nminvar = 2    
    acset = filter_acset(acset, nmincells, nminvar)
    lapply(acset, dim) #27,410, 163
    length(unique(acset[['featdata']][, 'feat']))
    ##nmincells:3 -> 4,188 genes.
    ##nmincells:5 -> 3,155 genes.

    ##Dump
    ##saveRDS(acset, file = )
    
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
    res = list(phased_list = phased_list, params = paramset)
    ##status: fin (15.54 -> ~16.04, 80 cores)
    
    ##Dump
    saveRDS(res, file = phased_rds)
    
    res = readRDS(phased_rds)

    lapply(res[['phased_list']], '[[', 'time')
    res[['params']]
    j.acset = res[['phased_list']][[5]]
    lapply(j.acset, dim) #27,410, 163
    length(unique(j.acset[['featdata']][, 'feat'])) #3155

    
}
