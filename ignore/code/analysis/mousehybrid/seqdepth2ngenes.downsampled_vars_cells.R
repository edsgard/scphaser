
###Assess complexity under varying parameters
###NB: The paths are relative to the root of the scphaser git repo
###
###SYNOPSIS
###library('devtools')
###devtools::load_all()
###source("ignore/code/analysis/mousehybrid/seqdepth2ngenes.downsampled_vars_cells.R")


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
ac.dir = file.path(cloud.dir, 'projects/scphaser/nogit/data/mousehybrid/all_snps/downsampled')
prefilt_rds_dir = '../nogit/data/mousehybrid/all_snps/downsampled'

##OUT
out_rds_dir = '../nogit/data/mousehybrid/all_snps/downsampled/vars_cells'
out_pdf_dir = './ignore/res/mousehybrid/all_snps/downsampled/vars_cells/perf/pdf'


##*###
##Files
##*###

##IN
varsampled_acset_rds = file.path(prefilt_rds_dir, 'acset_filt.rds')

##OUT
filt_acset_rds = file.path(out_rds_dir, 'acset_filt.rds')
param2ngenes_rds = file.path(out_rds_dir, 'seqdepth2ngenes.replace_no.rds')

##Libs
source('./ignore/code/analysis/performance.R')
library('BiocParallel')
library('dplyr')
library('tidyr')

main <- function(){

    ##read data
    acset = readRDS(varsampled_acset_rds)
        

    ##*###
    ##Downsample cells
    ##*###
    n.cells = 163 #same as in Borel et al

    ##Downsample cells
    cells = colnames(acset[['altcount']])
    cells.sample = sample(cells, n.cells)
    acset = subset_cols(acset, cells.sample)
    
    ##Filter vars with 0 counts
    acset = filter_zerorow(acset)
    lapply(acset, dim) #8,868
    
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
    lapply(acset, dim) #8,422 x 336

    ##Filter variants on n.cells monoallelic and feats with < 2 j.vars
    nmincells = 5
    nminvar = 2    
    acset = filter_acset(acset, nmincells, nminvar)
    lapply(acset, dim) #3,619 x 336
    length(unique(acset[['featdata']][, 'feat'])) #1,310

    ##Dump
    dir.create(dirname(filt_acset_rds), recursive = TRUE)
    saveRDS(acset, file = filt_acset_rds)
    
    
    ##*###
    ##Randomly sample n counts
    ##*###
    acset = readRDS(filt_acset_rds)    

    altcount = acset[['altcount']]
    refcount = acset[['refcount']]
    totcount = altcount + refcount

    ##number of reads
    ntot = sum(totcount) #3,280,983
    ncount = c(seq(2.5e5, ntot, 2.5e5), ntot)
    
    ##permutation iterations
    npermiter = 10 #10
    perm_iter = 1:npermiter
    
    ##specify paramset as all possible combinations of the params
    paramset = expand.grid(perm_iter, ncount, stringsAsFactors = FALSE)
    colnames(paramset) = c('perm_iter', 'ncount')
    
    ##parallelization
    ncores = 80
    bp_param = BiocParallel::MulticoreParam(workers = ncores)

    ##Filter vars and feats
    nparamset = nrow(paramset)
    acset_list = BiocParallel::bplapply(1:nparamset, filter_acset_par, BPPARAM = bp_param, paramset = paramset, acset = acset)
    ##status: sub (10.53 -> 10.54, 80 cores)

    ##get number of genes
    ngenes = unlist(lapply(acset_list, function(j.acset){length(unique(j.acset[['featdata']][, 'feat']))}))
    param2ngenes = cbind(paramset, ngenes)

    ##get number of variants
    nvars = unlist(lapply(acset_list, function(j.acset){nrow(j.acset[['featdata']])}))
    param2ngenes = cbind(param2ngenes, nvars)

    ##get fraction of genes
    n.bg.genes = 20268 ##see log.sh
    frac.genes = param2ngenes[, 'ngenes'] / n.bg.genes
    param2ngenes = cbind(param2ngenes, frac.genes)

    ##Dump
    saveRDS(param2ngenes, file = param2ngenes_rds)
    ##status: fin

    summary(param2ngenes[, 'ngenes'])

    
    ##*###
    ##Plot
    ##*###
    param2ngenes = readRDS(param2ngenes_rds)

    library('ggplot2')
    
    gg = ggplot(param2ngenes, aes_string(x = 'ncount', y = 'ngenes'))
    gg = gg + geom_line(stat = 'summary', fun.y = 'mean')
    gg = gg + geom_errorbar(stat = 'summary', fun.data = mean_se) #width = errbar.w

    ##tick breaks
    ##gg = gg + coord_cartesian(ylim = c(800, 3050))
    ##gg = gg + scale_y_continuous(breaks = seq(1000, 3000, 500))
    gg = gg + scale_x_continuous(breaks = c(seq(50, 336, 50), 336))
    
    ##Background
    gg = gg + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(panel.background = element_blank())
    gg = gg + theme(axis.text = element_text(colour="black"), axis.ticks = element_line(colour = 'black'))
    gg = gg + xlab('Number of reads')
    gg = gg + ylab('Number of phasable genes')

    j.pdf = file.path(out_pdf_dir, 'ncount2ngenes.replace_no.pdf')
    dir.create(dirname(j.pdf), recursive = TRUE)
    pdf(j.pdf)
    plot(gg)
    dev.off()

    
    gg = ggplot(param2ngenes, aes_string(x = 'ncount', y = 'frac.genes'))
    gg = gg + geom_line(stat = 'summary', fun.y = 'mean')
    gg = gg + geom_errorbar(stat = 'summary', fun.data = mean_se) #width = errbar.w

    ##tick breaks
    ##gg = gg + coord_cartesian(ylim = c(800, 3050))
    ##gg = gg + scale_y_continuous(breaks = seq(1000, 3000, 500))
    gg = gg + scale_x_continuous(breaks = c(seq(50, 336, 50), 336))
    
    ##Background
    gg = gg + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(panel.background = element_blank())
    gg = gg + theme(axis.text = element_text(colour="black"), axis.ticks = element_line(colour = 'black'))
    gg = gg + xlab('Number of reads')
    gg = gg + ylab('Number of phasable genes')
    gg = gg + theme(axis.line.x = element_line(color="black", size = 0.5), axis.line.y = element_line(color="black", size = 0.5))
    
    j.pdf = file.path(out_pdf_dir, 'ncount2fracgenes.replace_no.pdf')
    dir.create(dirname(j.pdf), recursive = TRUE)
    pdf(j.pdf)
    plot(gg)
    dev.off()

}

filter_acset_par <- function(j.param, paramset, acset, nmincells = 5, nminvar = 2){

    j_ncount = paramset[j.param, 'ncount']
    
    altcount = acset[['altcount']]
    refcount = acset[['refcount']]
    totcount = altcount + refcount
    ncount = sum(totcount)

    frac.down = j_ncount / ncount

    ##keep relative frequency btw alt and ref (skip sampling in this step)
    alt.frac = sum(altcount) / ncount
    ref.frac = 1 - alt.frac
    alt.frac.down = alt.frac * frac.down
    ref.frac.down = ref.frac * frac.down
    
    ##downsample
    altcount.down = downsample.matrix(altcount, alt.frac.down)
    refcount.down = downsample.matrix(refcount, ref.frac.down)
    
    ##set to downsampled data
    acset[['refcount']] = refcount.down
    acset[['altcount']] = altcount.down

    ##Call gt
    min_acount = 3
    fc = 3 #75/25
    acset = call_gt(acset, min_acount, fc)
    
    ##filter vars and genes
    acset_filt = filter_acset(acset, nmincells, nminvar)

    return(acset_filt)
}

downsample.matrix <- function(count.mat, frac.down, pseudo.count = 0){

    ##keep relative frequency between cols (skip sampling in this step)
    col2ncount = apply(count.mat, 2, sum)
    col2ncount.down = ceiling(col2ncount * frac.down)

    ##add pseudocount (skip, since we are down-sampling there is no need to do this)
    count.mat = count.mat + pseudo.count
    col2ncount = apply(count.mat, 2, sum)

    ##proportion of each element in a col
    probs.bycol = t(count.mat) / col2ncount
    ##each row: probs for a col

    ##for each col downsample    
    ##count.down = apply(Hmisc::rMultinom(probs.bycol, col2ncount.down), 1, table) ##Hmisc::rMultinom does not accept different sizes ("m = col2ncount.down").
    ##rownames(count.down) = rownames(count.mat)
    ##colnames(count.down) = colnames(count.mat)
    
    count.down = matrix(nrow = nrow(count.mat), ncol = ncol(count.mat), dimnames = list(rownames(count.mat), colnames(count.mat)))
    n.cols = length(col2ncount)
    for(j.col in 1:n.cols){
        jcol.ncount.down = col2ncount.down[j.col]
        jcol.probs = probs.bycol[j.col, ]
        jcol.count.down = rmultinom(1, jcol.ncount.down, prob = jcol.probs)
        count.down[, j.col] = jcol.count.down
    }

    return(count.down)
}
