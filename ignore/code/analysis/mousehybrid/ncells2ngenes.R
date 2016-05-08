
###Assess complexity under varying parameters
###NB: The paths are relative to the root of the scphaser git repo
###
###SYNOPSIS
###library('devtools')
###devtools::load_all()
###source("ignore/code/analysis/mousehybrid/ncells2ngenes.R")


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
out_rds_dir = '../nogit/data/mousehybrid'
out_pdf_dir = './ignore/res/mousehybrid/perf/pdf'


##*###
##Files
##*###

##IN
ref.counts.rds = file.path(cloud.dir, 'projects/scphaser/nogit/data/mousehybrid', paste('ref', '.counts.rds', sep = ''))
alt.counts.rds = file.path(cloud.dir, 'projects/scphaser/nogit/data/mousehybrid', paste('alt', '.counts.rds', sep = ''))
snp.annot.rds = file.path(cloud.dir, 'projects/scphaser/nogit/data/mousehybrid', 'snp.annot.rds')

##OUT
filt_acset_rds = file.path(out_rds_dir, 'acset_filt.rds')
ncells2ngenes_rds = file.path(out_rds_dir, 'ncells2ngenes.rds')

##Libs
source('./ignore/code/analysis/performance.R')
library('BiocParallel')
library('dplyr')
library('tidyr')

main <- function(){


    ##Read data
    load('../nogit/data/mousehybrid/mousehybrid.rda')

    ##create acset
    acset = new_acset(mousehybrid$featdata, mousehybrid$refcount, mousehybrid$altcount, mousehybrid$phenodata)
    lapply(acset, dim) #167,443 x 336
    length(unique(acset[['featdata']][, 'feat'])) #11,590
    feat2nvars = table(acset[['featdata']][, 'feat'])
    table(feat2nvars) #

    ##Filter vars with 0 counts
    acset = filter_zerorow(acset)
    lapply(acset, dim) #167,443
    
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
    lapply(acset, dim) #104,933 x 336

    ##Filter variants on n.cells monoallelic and feats with < 2 j.vars
    nmincells = 5
    nminvar = 2    
    acset = filter_acset(acset, nmincells, nminvar)
    lapply(acset, dim) #78,253 x 336
    length(unique(acset[['featdata']][, 'feat'])) #8,563

    ##Dump
    saveRDS(acset, file = filt_acset_rds)
    
    
    ##*###
    ##Randomly sample n cells
    ##*###
    acset = readRDS(filt_acset_rds)
    
    ##permutation iterations
    npermiter = 10 #10
    perm_iter = 1:npermiter

    ##number of cells
    ncells = c(seq(25, 336, 25), 336)    
    
    ##specify paramset as all possible combinations of the params
    paramset = expand.grid(perm_iter, ncells, stringsAsFactors = FALSE)
    colnames(paramset) = c('perm_iter', 'ncells')
    
    ##parallelization
    ncores = 80
    bp_param = BiocParallel::MulticoreParam(workers = ncores)

    ##Filter vars and feats
    nparamset = nrow(paramset)
    acset_list = BiocParallel::bplapply(1:nparamset, filter_acset_par, BPPARAM = bp_param, paramset = paramset, acset = acset)
    ##status: sub (11.47 -> , 80 cores)

    ##get number of genes
    ngenes = unlist(lapply(acset_list, function(j.acset){length(unique(j.acset[['featdata']][, 'feat']))}))
    ncells2ngenes = cbind(paramset, ngenes)

    ##get number of variants
    nvars = unlist(lapply(acset_list, function(j.acset){nrow(j.acset[['featdata']])}))
    ncells2ngenes = cbind(ncells2ngenes, nvars)
    
    ##Dump
    saveRDS(ncells2ngenes, file = ncells2ngenes_rds)

    
    ##*###
    ##Plot
    ##*###
    ncells2ngenes = readRDS(ncells2ngenes_rds)

    library('ggplot2')
    
    gg = ggplot(ncells2ngenes, aes_string(x = 'ncells', y = 'ngenes'))
    gg = gg + geom_line(stat = 'summary', fun.y = 'mean')
    gg = gg + geom_errorbar(stat = 'summary', fun.data = mean_se) #width = errbar.w

    ##tick breaks
    ##gg = gg + coord_cartesian(ylim = c(800, 3050))
    ##gg = gg + scale_y_continuous(breaks = seq(1000, 3000, 500))
    gg = gg + scale_x_continuous(breaks = c(seq(50, 336, 50), 336))
    
    ##Background
    gg = gg + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(panel.background = element_blank())
    gg = gg + theme(axis.text = element_text(colour="black"), axis.ticks = element_line(colour = 'black'))
    gg = gg + xlab('Number of sequenced cells')
    gg = gg + ylab('Number of phasable genes')

    j.pdf = file.path(out_pdf_dir, 'ncells2ngenes.resampled.pdf')
    dir.create(dirname(j.pdf), recursive = TRUE)
    pdf(j.pdf)
    plot(gg)
    dev.off()
    
}

filter_acset_par <- function(j.param, paramset, acset, nmincells = 5, nminvar = 2){

    j_ncells = paramset[j.param, 'ncells']
    
    ##subset cells
    phenodata = acset[['phenodata']]
    samples = phenodata[, 'sample']
    j_samples = sample(samples, j_ncells, replace = TRUE)
    acset_filt = subset_cols(acset, j_samples)

    ##filter vars and genes
    acset_filt = filter_acset(acset_filt, nmincells, nminvar)

    return(acset_filt)
}
