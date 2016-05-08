
###Assess complexity under varying parameters
###NB: The paths are relative to the root of the scphaser git repo
###
###SYNOPSIS
###library('devtools')
###devtools::load_all()
###source("ignore/code/analysis/borel/ase2nvars.R")


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
out_pdf_dir = './ignore/res/borel/perf/pdf'


##*###
##Files
##*###

##IN
ref.counts.rds = file.path(cloud.dir, 'projects/scphaser/nogit/data/borel', paste('ref', '.counts.rds', sep = ''))
alt.counts.rds = file.path(cloud.dir, 'projects/scphaser/nogit/data/borel', paste('alt', '.counts.rds', sep = ''))
snp.annot.rds = file.path(cloud.dir, 'projects/scphaser/nogit/data/borel', 'snp.annot.rds')

##OUT
filt_acset_rds = file.path(out_rds_dir, 'acset_filt_homo.rds')

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


    ##exclude entries with <= 2 reads
    mincount = 3

    totcount = acset[['altcount']] + acset[['refcount']]
    ase = acset[['altcount']] / totcount
    cells = colnames(ase)
    ase.list = list()
    for(j.cell in cells){
        pass.var = rownames(totcount)[which(totcount[, j.cell] >= mincount)]
        pass.ase = ase[pass.var, j.cell]
        ase.list[[j.cell]] = pass.ase
    }
    
    ##For each cell calculate the ase2nvars density
    ase2nvar.dens.list = lapply(ase.list, density, from = 0, to = 1)

    ##double-check that breaks are identical
    x.list = lapply(ase2nvar.dens.list, '[[', 'x')
    a = do.call(rbind, x.list)
    b = lapply(as.data.frame(a), unique)
    unique(unlist(lapply(b, length)))

    ##list -> long
    ##ase, nvar.dens, cell
    cells = names(ase2nvar.dens.list)
    ase2nvar.dens.df = matrix(nrow = 0, ncol = 3)
    for(j.cell in cells){
        j.dens = ase2nvar.dens.list[[j.cell]]
        ase = j.dens[['x']]
        nvar.dens = j.dens[['y']]
        n.breaks = length(ase)
        cell = rep(j.cell, n.breaks)
        j.mat = cbind(ase, nvar.dens, cell)
        ase2nvar.dens.df = rbind(ase2nvar.dens.df, j.mat)
    }
    ase2nvar.dens.df = as.data.frame(ase2nvar.dens.df, stringsAsFactors = FALSE)
    ase2nvar.dens.df[, 'ase'] = as.numeric(ase2nvar.dens.df[, 'ase'])
    ase2nvar.dens.df[, 'nvar.dens'] = as.numeric(ase2nvar.dens.df[, 'nvar.dens'])

    ##Dump
    saveRDS(ase2nvar.dens.df, file = file.path(out_rds_dir, 'ase2nvars.dens.rds'))

    
    ##*###
    ##Plot
    ##*###
    library('ggplot2')
    
    gg = ggplot(ase2nvar.dens.df, aes_string(x = 'ase', y = 'nvar.dens'))
    gg = gg + geom_line(stat = 'summary', fun.y = 'median')
    ##gg = gg + stat_summary(geom = 'ribbon', fun.data = mean_se)
    ##gg = gg + stat_summary(geom = 'ribbon', fun.ymin = function(x){mean(x) - sd(x)}, fun.ymax = function(x){mean(x) + sd(x)})
    gg = gg + stat_summary(geom = 'ribbon', fun.ymin = function(x){quantile(x, probs = 0.25)}, fun.ymax = function(x){quantile(x, probs = 0.75)})
    
    ##Background
    gg = gg + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(panel.background = element_blank())
    gg = gg + theme(axis.text = element_text(colour="black"), axis.ticks = element_line(colour = 'black'))
    gg = gg + xlab('Allele specific expression')
    gg = gg + ylab('Number of variants (density)')

    plot(gg)
    
    
}
