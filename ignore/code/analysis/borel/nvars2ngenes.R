
###Assess complexity under varying parameters
###NB: The paths are relative to the root of the scphaser git repo
###
###SYNOPSIS
###library('devtools')
###devtools::load_all()
###source("ignore/code/analysis/borel/nvars2ngenes.R")


##Params
sys = 'rs13'
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

ac.dirs = list.dirs(file.path(cloud.dir, sprintf('projects/scphaser/nogit/data/borel')), recursive = FALSE)
ac.dirs = c(ac.dirs, dirname(ac.dirs[1]))

##OUT
out_rds_dir = '../nogit/data/borel'
out_pdf_dir = './ignore/res/borel/pdf'


##*###
##Files
##*###

##OUT
##acset_list_rds = file.path(out_rds_dir, 'acset_list.rds')
param2ngenes_rds = file.path(out_rds_dir, 'nvars2ngenes.rds')

##Libs
source('./ignore/code/analysis/performance.R')
library('BiocParallel')
library('dplyr')
library('tidyr')

main <- function(){

    acset_list = list()
    for(j.dir in ac.dirs){

        print(j.dir)
        
        ref.counts.rds = file.path(j.dir, 'ref.counts.rds')
        alt.counts.rds = file.path(j.dir, 'alt.counts.rds')
        snp.annot.rds = file.path(j.dir, 'snp.annot.rds')
        
        ##Read data
        altcount = readRDS(alt.counts.rds)
        refcount = readRDS(ref.counts.rds)
        featdata = readRDS(snp.annot.rds)
        
        ##create acset
        acset = new_acset(featdata, refcount, altcount)
        lapply(acset, dim) #398,399 x 336
        length(unique(acset[['featdata']][, 'feat'])) #17,790
        
        ##Filter vars with 0 counts
        acset = filter_zerorow(acset)
        lapply(acset, dim) #167,443
        length(unique(acset[['featdata']][, 'feat']))

        ##Call gt
        min_acount = 3
        fc = 3 #75/25
        acset = call_gt(acset, min_acount, fc)

        ##Filter variants on number of cells where ASE towards the same allele (as to remove probable false positive calls)
        alpha = 0.1
        mono.ase = 0.1
        if(!(mono.ase == 0)){
            acset = filter_homovars(acset, alpha = alpha, mono.ase = mono.ase)
        }
        lapply(acset, dim) #104,933 x 336
        
        ##filter vars and genes
        nmincells = 5
        nminvar = 2
        acset_filt = filter_acset(acset, nmincells, nminvar)

        acset_list[[j.dir]] = acset_filt
    }
    names(acset_list) = basename(names(acset_list))
    names(acset_list) = sub('borel', 'N787883', names(acset_list))
    
    ##Dump
    ##saveRDS(acset_list, file = acset_list_rds)
    
    
    ##get number of genes
    ngenes = unlist(lapply(acset_list, function(j.acset){length(unique(j.acset[['featdata']][, 'feat']))}))
    nvars.in = as.integer(sub('.*N', '', names(acset_list)))
    param2ngenes = cbind(ngenes, nvars.in)

    ##get number of variants
    nvars = unlist(lapply(acset_list, function(j.acset){nrow(j.acset[['featdata']])}))
    param2ngenes = cbind(param2ngenes, nvars)

    ##get fraction of genes
    n.bg.genes = 15556 ##see ngenes.R
    frac.genes = param2ngenes[, 'ngenes'] / n.bg.genes
    param2ngenes = cbind(param2ngenes, frac.genes)
    param2ngenes = as.data.frame(param2ngenes)
    
    ##Dump
    saveRDS(param2ngenes, file = param2ngenes_rds)
    ##status: fin

    summary(param2ngenes[, 'ngenes'])

    
    ##*###
    ##Plot
    ##*###
    param2ngenes = readRDS(param2ngenes_rds)

    library('ggplot2')

    
    gg = ggplot(param2ngenes, aes_string(x = 'nvars.in', y = 'ngenes'))
    gg = gg + geom_line(stat = 'summary', fun.y = 'mean')

    ##tick breaks
    ##gg = gg + coord_cartesian(ylim = c(800, 3050))
    ##gg = gg + scale_y_continuous(breaks = seq(1000, 3000, 500))
    ##gg = gg + scale_x_continuous(breaks = c(seq(25e6, 125e6, 25e6)))
    
    ##Background
    gg = gg + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(panel.background = element_blank())
    gg = gg + theme(axis.text = element_text(colour="black", size = 15), axis.ticks = element_line(colour = 'black'), axis.title = element_text(size = 15))
    gg = gg + xlab('Number of heterozygous variants')
    gg = gg + ylab('Number of phasable genes')
    gg = gg + theme(axis.line.x = element_line(color="black", size = 0.5), axis.line.y = element_line(color="black", size = 0.5))

    pdf.h = 5
    pdf.w = 5
    
    j.pdf = file.path(out_pdf_dir, 'nvars2ngenes.pdf')
    dir.create(dirname(j.pdf), recursive = TRUE)
    pdf(j.pdf, height = pdf.h, width = pdf.w)
    plot(gg)
    dev.off()


    ##*###
    ##y: frac.genes
    ##*###
    gg = ggplot(param2ngenes, aes_string(x = 'nvars.in', y = 'frac.genes'))
    gg = gg + geom_line(stat = 'summary', fun.y = 'mean')

    ##tick breaks
    ##gg = gg + coord_cartesian(ylim = c(800, 3050))
    ##gg = gg + scale_y_continuous(breaks = seq(1000, 3000, 500))
    ##gg = gg + scale_x_continuous(breaks = c(seq(50, 336, 50), 336))
    
    ##Background
    gg = gg + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(panel.background = element_blank())
    gg = gg + theme(axis.text = element_text(colour="black", size = 15), axis.ticks = element_line(colour = 'black'), axis.title = element_text(size = 15))
    gg = gg + xlab('Number of heterozygous variants')
    gg = gg + ylab('Number of phasable genes')
    gg = gg + theme(axis.line.x = element_line(color="black", size = 0.5), axis.line.y = element_line(color="black", size = 0.5))

    pdf.h = 5
    pdf.w = 5
 
    j.pdf = file.path(out_pdf_dir, 'nvars2fracgenes.pdf')
    dir.create(dirname(j.pdf), recursive = TRUE)
    pdf(j.pdf, height = pdf.h, width = pdf.w)
    plot(gg)
    dev.off()
}
